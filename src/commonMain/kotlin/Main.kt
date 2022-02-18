import space.kscience.kmath.linear.DoubleLinearSpace.dot
import space.kscience.kmath.nd.MutableStructure2D
import space.kscience.kmath.nd.Structure2D
import space.kscience.kmath.nd.as1D
import space.kscience.kmath.nd.as2D
import space.kscience.kmath.tensors.api.Tensor
import space.kscience.kmath.tensors.core.*
import kotlin.math.abs
import kotlin.math.min

object Main {

    fun runJacobi() {
        val n = 2
        val src = DoubleTensorAlgebra.fromArray(intArrayOf(n, n), doubleArrayOf(1.0, 2.0, 2.0, 4.0))
        println(src)
        val (values, vectors) = src.symEigJacobi(DoubleTensorAlgebra, 0.001)
        println(values)
        println(vectors)
        DoubleTensorAlgebra.withBroadcast {
            val sv = (src.dot(vectors.transpose()))
            val vv = (values * vectors.transpose())
            println(sv.eq(vv, 0.001))
        }
    }

    fun DoubleTensorAlgebra.checkSymmetric(
        tensor: Tensor<Double>, epsilon: Double = 1e-6
    ) =
        check(tensor.eq(tensor.transpose(), epsilon)) {
            "Tensor is not symmetric about the last 2 dimensions at precision $epsilon"
        }

    fun Tensor<Double>.symEigJacobi(algebra: DoubleTensorAlgebra, epsilon: Double): Pair<DoubleTensor, DoubleTensor> {
        with(algebra) {
            check(shape.size == 2)
            checkSymmetric(this@symEigJacobi, epsilon)
            val S = copy().as2D()
            val n = S.shape[0]
            fun indexOfMaxInRow(row: Int): Int {
                var res = if (row > 0) 0 else 1
                for (i in 0 until n) {
                    if (i != row && abs(S[row, i]) > abs(S[row, res])) {
                        res = i
                    }
                }
                return res
            }

            val eigenvalues = fromArray(intArrayOf(n), DoubleArray(n) { S[it, it] })
            val eigenvectors = eye(n)

            if (n <= 1) {
                return eigenvalues to eigenvectors
            }

            val eigenvaluesW = eigenvalues.as1D()
            val eigenvectorsW = eigenvectors.as2D()

            val changed = BooleanArray(n) { true }
            var changedCount = n

            val maxind = IntArray(n) { indexOfMaxInRow(it) }


            fun updateEigenvalues(k: Int, t: Double) {
                eigenvaluesW[k] += t
                when {
                    changed[k] && 2 * n * abs(t) <= epsilon -> {
                        changed[k] = false
                        changedCount--
                    }
                    !changed[k] && 2 * n * abs(t) > epsilon -> {
                        changed[k] = true
                        changedCount++
                    }
                }
            }

            fun pseudoSign(x: Double) = if (x < 0) -1 else 1

            while (changedCount > 0) {
                var pivotRow = 0
                for (i in 1 until n) {
                    if (abs(S[i, maxind[i]]) > abs(S[pivotRow, maxind[pivotRow]])) {
                        pivotRow = i
                    }
                }

                val pivotCol = maxind[pivotRow]
                val pivot = S[pivotRow, pivotCol]

                val y = (eigenvaluesW[pivotCol] - eigenvaluesW[pivotRow]) / 2
                val d = abs(y) + kotlin.math.sqrt(pivot * pivot + y * y)
                val r = kotlin.math.sqrt(pivot * pivot + d * d)

                val c = d / r
                val s = pseudoSign(y) * pivot / r
                val t = pseudoSign(y) * (pivot * pivot) / d

                updateEigenvalues(pivotCol, t)
                updateEigenvalues(pivotRow, -t)

                fun MutableStructure2D<Double>.rotate(k: Int, l: Int, i: Int, j: Int) {
                    val skl = this[k, l]
                    val sij = this[i, j]
                    this[k, l] = c * skl - s * sij
                    this[i, j] = s * skl + c * sij
                }

                for (i in 0 until n) {
                    S.rotate(i, pivotRow, i, pivotCol)
                }
                for (i in 0 until n) {
                    S.rotate(pivotRow, i, pivotCol, i)
                    eigenvectorsW.rotate(pivotRow, i, pivotCol, i)
                }
                maxind[pivotRow] = indexOfMaxInRow(pivotRow)
                maxind[pivotCol] = indexOfMaxInRow(pivotCol)
            }
            return eigenvalues to eigenvectors
        }
    }


    fun runStrassen() {
        val a = DoubleTensorAlgebra.fromArray(
            intArrayOf(3, 3), doubleArrayOf(
                1.0, 2.0, 3.0,
                5.0, 6.0, 7.0,
                1.0, 2.0, 3.0
            )
        )
        val b = DoubleTensorAlgebra.fromArray(
            intArrayOf(3, 3), doubleArrayOf(
                1.0, 2.0, 3.0,
                5.0, 6.0, 7.0,
                1.0, 2.0, 3.0
            )
        )
        DoubleTensorAlgebra.withBroadcast {
            val d = (a dot b).as2D()
            val sd = (a.as2D() strassenDot b.as2D())
            check(d.eq(sd))
        }
    }

    private fun Tensor<Double>.strassenValid(): Boolean {
        return (shape.size == 2 &&
                shape[0] == shape[1])
    }

    infix fun Tensor<Double>.strassenDot(other: Tensor<Double>): Tensor<Double> {
        check(this.strassenValid())
        check(other.strassenValid())
        check(this.shape[1] == other.shape[0])
        return strassenDoInternal(this, other)
    }

    private fun strassenDoInternal(a: Tensor<Double>, b: Tensor<Double>): Tensor<Double> =
        DoubleTensorAlgebra.withBroadcast {
            val n = a.shape[0]
            if (n <= 2) {
                return@withBroadcast (a dot b)
            }
            val m = (n + 1) / 2
            val k = m * 2
            val a11 = a.submatrix(0, 0, m, m)
            val a12 = a.submatrix(0, m, m, k)
            val a21 = a.submatrix(m, 0, k, m)
            val a22 = a.submatrix(m, m, k, k)

            val b11 = b.submatrix(0, 0, m, m)
            val b12 = b.submatrix(0, m, m, k)
            val b21 = b.submatrix(m, 0, k, m)
            val b22 = b.submatrix(m, m, k, k)

            val m1 = strassenDoInternal(a11 + a22, b11 + b22)
            val m2 = strassenDoInternal(a21 + a22, b11)
            val m3 = strassenDoInternal(a11, b12 - b22)
            val m4 = strassenDoInternal(a22, b21 - b11)
            val m5 = strassenDoInternal(a11 + a12, b22)
            val m6 = strassenDoInternal(a21 - a11, b11 + b12)
            val m7 = strassenDoInternal(a12 - a22, b21 + b22)

            DoubleTensorAlgebra.zero(n, n).also { c ->
                put(c, m1 + m4 - m5 + m7, 0, 0, m, m)
                put(c, m3 + m5, 0, m, m, k)
                put(c, m2 + m4, m, 0, k, m)
                put(c, m1 - m2 + m3 + m6, m, m, k, k)
            }
        }

    fun put(m: Tensor<Double>, x: Tensor<Double>, x1: Int, y1: Int, x2: Int, y2: Int) {
        for (i in 0 until min(x2 - x1, m.shape[0] - x1)) {
            for (j in 0 until min(y2 - y1, m.shape[1] - y1)) {
                m[intArrayOf(x1 + i, y1 + j)] = x[intArrayOf(i, j)]
            }
        }
    }

    fun Tensor<Double>.submatrix(a1: Int, b1: Int, a2: Int, b2: Int): MatrixView {
        return if (this is MatrixView) {
            MatrixView(data, x1 + a1, y1 + b1, x1 + a2, y1 + b2)
        } else {
            MatrixView(this, a1, b1, a2, b2)
        }
    }

    class MatrixView(
        val data: Tensor<Double>,
        val x1: Int,
        val y1: Int,
        val x2: Int,
        val y2: Int
    ) : Tensor<Double> {

        operator fun set(i: Int, j: Int, value: Double) {
            val x = x1 + i
            val y = y1 + j
            if (x < data.shape[0] && y < data.shape[1]) {
                data[intArrayOf(x, y)] = value
            }
        }

        override fun set(index: IntArray, value: Double) {
            this[index[0], index[1]] = value
        }

        operator fun get(i: Int, j: Int): Double {
            val x = x1 + i
            val y = y1 + j
            if (x >= data.shape[0] || y >= data.shape[1]) {
                return 0.0
            }
            return data[intArrayOf(x, y)]
        }

        override val shape by lazy { intArrayOf(x2 - x1, y2 - y1) }

        override fun elements() = sequence {
            for (i in 0 until shape[0])
                for (j in 0 until shape[1]) yield(intArrayOf(i, j) to get(i, j))
        }

        override fun get(index: IntArray): Double {
            return this[index[0], index[1]]
        }

        override fun toString(): String {
            return (0 until shape[0]).joinToString(", ", "[", "]") { i ->
                (0 until shape[1]).joinToString(", ", "[", "]") { j ->
                    get(i, j).toString()
                }
            }
        }
    }
}