import Main.symEigJacobi
import space.kscience.kmath.tensors.api.Tensor
import space.kscience.kmath.tensors.core.DoubleTensorAlgebra
import space.kscience.kmath.tensors.core.withBroadcast
import kotlin.random.Random
import kotlin.random.nextInt
import kotlin.test.Test
import kotlin.test.assertTrue

class JacobiTest {

    private fun DoubleTensorAlgebra.validate(m: Tensor<Double>, epsilon: Double): Boolean = withBroadcast {
        val (values, vectors) = m.symEigJacobi(this, epsilon)
        val sv = (m.dot(vectors.transpose()))
        val vv = (values * vectors.transpose())
        return@withBroadcast sv.eq(vv, values.max() * epsilon)
    }

    @Test
    fun testJacobi() {
        val epsilon = 0.001
        val maxSize = 50
        val elemScale = 1000
        with(DoubleTensorAlgebra) {
            val rng = Random(6741)
            repeat(1000) { iter ->
                val n = rng.nextInt(1..maxSize)
                val m = zeros(intArrayOf(n, n))
                for (i in 0 until n) {
                    for (j in i until n) {
                        m[intArrayOf(i, j)] = rng.nextDouble() * elemScale
                        m[intArrayOf(j, i)] = m[intArrayOf(i, j)]
                    }
                }
                assertTrue(validate(m, epsilon), "Error on iter=$iter")
            }
        }
    }

}