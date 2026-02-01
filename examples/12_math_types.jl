# Math Types Example
# Demonstrates Polynom, LinearSpline, and CubicSpline from Phase 10.
# These are used internally for density profiles and interpolation.

using PROPOSAL

println("=== Math Types ===\n")

# --- Polynom ---
println("--- Polynom ---")
# Coefficients: p(x) = 1 + 2x + 3x^2
p = Polynom([1.0, 2.0, 3.0])
println("p(x) = 1 + 2x + 3x^2")
for x in [0.0, 0.5, 1.0, 2.0, -1.0]
    expected = 1.0 + 2.0*x + 3.0*x^2
    actual = evaluate(p, x)
    println("  p($x) = $actual (expected: $expected)")
end

# Derivative: p'(x) = 2 + 6x
dp = derive(p)
println("\np'(x) = 2 + 6x")
for x in [0.0, 1.0, 2.0]
    println("  p'($x) = $(evaluate(dp, x))")
end

# Antiderivative: P(x) = C + x + x^2 + x^3
ap = antiderivative(p, 0.0)
println("\nP(x) = x + x^2 + x^3 (with C=0)")
for x in [0.0, 1.0, 2.0]
    println("  P($x) = $(evaluate(ap, x))")
end

# Get coefficients
coeffs = get_coefficient(p)
println("\nCoefficients: $coeffs")

# --- Linear Spline ---
println("\n--- Linear Spline ---")
xs = [0.0, 1.0, 2.0, 3.0, 4.0]
ys = [0.0, 1.0, 4.0, 9.0, 16.0]  # y = x^2 sampled
ls = LinearSpline(xs, ys)
println("Linear interpolation of y=x^2:")
for x in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    println("  f($x) = $(evaluate(ls, x)) (exact: $(x^2))")
end

# --- Cubic Spline ---
println("\n--- Cubic Spline ---")
cs = CubicSpline(xs, ys)
println("Cubic interpolation of y=x^2:")
for x in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    println("  f($x) = $(round(evaluate(cs, x), sigdigits=6)) (exact: $(x^2))")
end

# --- Linear spline of a linear function should be exact ---
println("\n--- Verification: linear spline of f(x)=x+1 ---")
xs2 = [0.0, 1.0, 2.0, 3.0]
ys2 = [1.0, 2.0, 3.0, 4.0]
ls2 = LinearSpline(xs2, ys2)
for x in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    val = evaluate(ls2, x)
    expected = x + 1.0
    err = abs(val - expected)
    println("  f($x) = $val (expected: $expected, err: $err)")
end

println("\nDone.")
