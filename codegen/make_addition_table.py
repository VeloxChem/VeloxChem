"""Builds the real solid harmonics symbolically with the recursion of
make_solid_harmonics.py, and solves exactly for the coefficients of the addition
theorem which splits S_l,m(a + b) by bidegree.

Set VLX_HARM_PROBE to a program which prints the harmonics the library computes,
and the symbolic recursion is checked against it before the table is built. See
the section on this generator in README.md."""
import os, subprocess, sys, random
import sympy as sp

x, y, z = sp.symbols('x y z')
ax, ay, az, bx, by, bz = sp.symbols('ax ay az bx by bz')

_cache = {}
def solid(l, m, X=x, Y=y, Z=z):
    """Real regular solid harmonic, by the recursion of make_solid_harmonics.py."""
    key = (l, m)
    if key not in _cache:
        if l == 0:
            e = sp.Integer(1)
        elif l == 1:
            e = {-1: y, 0: z, 1: x}[m]
        else:
            r2 = x * x + y * y + z * z
            if m == l:
                f = sp.sqrt(sp.Rational(2 * l - 1, 2 * l))
                e = f * (x * solid(l - 1, l - 1) - y * solid(l - 1, -(l - 1)))
            elif m == -l:
                f = sp.sqrt(sp.Rational(2 * l - 1, 2 * l))
                e = f * (y * solid(l - 1, l - 1) + x * solid(l - 1, -(l - 1)))
            elif abs(m) == l - 1:
                e = sp.sqrt(sp.Integer(2 * l - 1)) * z * solid(l - 1, m)
            else:
                fz = sp.sqrt(sp.Rational((2 * l - 1) ** 2, (l + m) * (l - m)))
                fr = sp.sqrt(sp.Rational((l - 1 + m) * (l - 1 - m), (l + m) * (l - m)))
                e = fz * z * solid(l - 1, m) - fr * r2 * solid(l - 2, m)
        _cache[key] = sp.expand(e)
    return _cache[key].subs({x: X, y: Y, z: Z}, simultaneous=True)

PROBE = os.environ.get('VLX_HARM_PROBE', '')

def check_convention(lmax=6, ntrials=4):
    if not PROBE:
        print('VLX_HARM_PROBE is not set, the convention check is skipped')
        return 0.0
    random.seed(11)
    worst = 0.0
    for _ in range(ntrials):
        v = [random.uniform(-2, 2) for _ in range(3)]
        out = subprocess.run([PROBE, str(lmax)] + ['%.17g' % t for t in v] + ['0', '0', '0'],
                             capture_output=True, text=True).stdout
        for line in out.strip().splitlines():
            l, m, val = line.split()
            mine = float(solid(int(l), int(m)).subs({x: v[0], y: v[1], z: v[2]}))
            worst = max(worst, abs(mine - float(val)))
    return worst

def table(l):
    """C[(l1, m1, m2)][m] for the bidegree split of S_{l,m}(a+b)."""
    out = {}
    for l1 in range(l + 1):
        l2 = l - l1
        prods = [((m1, m2), sp.expand(solid(l1, m1, ax, ay, az) * solid(l2, m2, bx, by, bz)))
                 for m1 in range(-l1, l1 + 1) for m2 in range(-l2, l2 + 1)]
        mons = sorted({mo for _, p in prods for mo in p.as_poly(ax, ay, az, bx, by, bz).monoms()})
        M = sp.Matrix([[p.as_poly(ax, ay, az, bx, by, bz).coeff_monomial(mo) for _, p in prods]
                       for mo in mons])
        for m in range(-l, l + 1):
            full = sp.expand(solid(l, m, ax + bx, ay + by, az + bz))
            # the part of bidegree (l1, l2): degree l1 in a, l2 in b
            part = sum(c * sp.prod([s ** e for s, e in zip((ax, ay, az, bx, by, bz), mo)])
                       for mo, c in full.as_poly(ax, ay, az, bx, by, bz).terms()
                       if sum(mo[:3]) == l1 and sum(mo[3:]) == l2)
            part = sp.expand(part)
            rhs = sp.Matrix([[part.as_poly(ax, ay, az, bx, by, bz).coeff_monomial(mo)] for mo in mons]) \
                  if part != 0 else sp.zeros(len(mons), 1)
            sol = M.solve_least_squares(rhs)
            assert sp.simplify(M * sol - rhs) == sp.zeros(len(mons), 1), 'no exact decomposition'
            for ((m1, m2), _), c in zip(prods, sol):
                c = sp.nsimplify(sp.simplify(c))
                if c != 0:
                    out.setdefault((l1, m1, m2), {})[m] = c
    return out

if __name__ == '__main__':
    w = check_convention()
    print('symbolic harmonics against the library, l <= 6 : max deviation %.2e' % w)
    assert w < 1.0e-12, 'the symbolic recursion is not the library convention'
    l = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    t = table(l)
    print('\nl = %d, non-zero coefficients:' % l)
    for (l1, m1, m2) in sorted(t):
        for m in sorted(t[(l1, m1, m2)]):
            print('  l1=%d m1=%+d  l2=%d m2=%+d  ->  m=%+d : %s'
                  % (l1, m1, l - l1, m2, m, t[(l1, m1, m2)][m]))
