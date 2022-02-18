import Main.strassenDot
import space.kscience.kmath.tensors.api.Tensor
import space.kscience.kmath.tensors.core.DoubleTensorAlgebra
import kotlin.math.exp
import kotlin.math.round
import kotlin.random.Random
import kotlin.random.nextInt
import kotlin.test.Test
import kotlin.test.assertTrue

class StrassenTest {

    fun DoubleTensorAlgebra.validate(a: Tensor<Double>, b: Tensor<Double>): Boolean {
        val d = a dot b
        val sd = a strassenDot b
        return d eq sd
    }

    private fun Random.nextTensor(shape: IntArray, elemScale: Double): Tensor<Double> =
        DoubleTensorAlgebra.produce(shape) { round(nextDouble() * elemScale) }

    @Test
    fun testStrassen() {
        val maxSize = 128
        val elemSize = 1000.0
        with(DoubleTensorAlgebra) {
            val rng = Random(4321)
            repeat(50) { iter ->
                val n = rng.nextInt(1..maxSize)

                val a = rng.nextTensor(intArrayOf(n, n), elemSize)
                val b = rng.nextTensor(intArrayOf(n, n), elemSize)
                assertTrue(validate(a, b), "Error on iter=$iter")
            }
        }
    }
}