plugins {
    kotlin("multiplatform") version "1.5.10"
}

group = "com.frgm"
version = "1.0-SNAPSHOT"
repositories {
    mavenCentral()
    maven("https://repo.kotlin.link")
}
dependencies {

}

kotlin {
    jvm {
        compilations.all {
            kotlinOptions.jvmTarget = "1.8"
        }
        testRuns["test"].executionTask.configure {
            useJUnit()
        }
    }
    js(IR) {
        browser {
            commonWebpackConfig {
                cssSupport.enabled = true
            }
        }
    }
    val hostOs = System.getProperty("os.name")
    val isMingwX64 = hostOs.startsWith("Windows")
    val nativeTarget = when {
        hostOs == "Mac OS X" -> null // no kmath for macos
        hostOs == "Linux" -> linuxX64("native")
        isMingwX64 -> mingwX64("native")
        else -> throw GradleException("Host OS is not supported in Kotlin/Native.")
    }

    sourceSets {
        val commonMain by getting
        val commonTest by getting
        val jvmMain by getting
        val jvmTest by getting
        val jsMain by getting
        val jsTest by getting
        val nativeMain by getting
        val nativeTest by getting

        commonMain.dependencies {
            api("space.kscience:kmath-core:0.3.0-dev-17")
            api("space.kscience:kmath-tensors:0.3.0-dev-17")
        }

        commonTest.dependencies {
            implementation(kotlin("test"))
        }
    }

}

tasks.withType<org.jetbrains.kotlin.gradle.tasks.KotlinCompile>().configureEach {
    kotlinOptions {
        freeCompilerArgs = listOf("-Xjvm-default=compatibility")
    }
}