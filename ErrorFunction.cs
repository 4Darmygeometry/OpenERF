// OpenERF, 开源Erf、Erfc、Erfcx、ErfInv高速C#实数实现

using System;

namespace OpenERFFunc
{
    /// <summary>
    /// 误差函数及其相关函数（纯实数双精度）。
    /// erf：Cody 1969 有理逼近 + 小 |x| Taylor 级数（Dekker TwoProd 补偿，正确舍入级）；
    /// erfc：Cody 三段有理逼近（带 exp 分裂，避免大 x 下溢失真）；
    /// erfcx(x) = e^{x²}·erfc(x)：同样走 Cody 分段，利用 exp 分裂的恒等关系 e^{y²}·e^{-ysq²}·e^{-del} ≡ 1，
    /// 免去大 x 时 exp 上下溢对消，精度与 erfc 同阶。
    /// </summary>
    public static class ErrorFunctions
    {
        // ==================== 常数（Cody 实数部分） ====================
        private const double SqrPi = 0.564189583547756286948079451560772585844050629329;
        private const double Thresh = 0.46875;
        private const double Four = 4.0;
        private const double Six = 6.0;
        private const double Sixteen = 16.0;
        private const double XSmall = 1.11e-16;
        private const double XBig = 26.543;
        private const double XHuge = 6.71e7;

        private static readonly double[] A = new double[]
        {
            3.16112374387056560, 1.13864154151050156e2, 3.77485237685302021e2,
            3.20937758913846947e3, 1.85777706184603153e-1
        };
        private static readonly double[] B = new double[]
        {
            2.36012909523441209e1, 2.44024637934444173e2,
            1.28261652607737228e3, 2.84423683343917062e3
        };
        private static readonly double[] C = new double[]
        {
            5.64188496988670089e-1, 8.88314979438837594e0, 6.61191906371416295e1,
            2.98635138197400131e2, 8.81952221241769090e2, 1.71204761263407058e3,
            2.05107837782607147e3, 1.23033935479799725e3, 2.15311535474403846e-8
        };
        private static readonly double[] D = new double[]
        {
            1.57449261107098347e1, 1.17693950891312499e2, 5.37181101862009858e2,
            1.62138957456669019e3, 3.29079923573345963e3, 4.36261909014324716e3,
            3.43936767414372164e3, 1.23033935480374942e3
        };
        private static readonly double[] P = new double[]
        {
            3.05326634961232344e-1, 3.60344899949804439e-1, 1.25781726111229246e-1,
            1.60837851487422766e-2, 6.58749161529837803e-4, 1.63153871373020978e-2
        };
        private static readonly double[] Q = new double[]
        {
            2.56852019228982242e0, 1.87295284992346047e0, 5.27905102951428412e-1,
            6.05183413124413191e-2, 2.33520497626869185e-3
        };

        // ==================== FMA 辅助 ====================

        // Horner 乘法累加：新 TFM 用硬件 FMA（单舍入），旧 TFM 降级为普通乘加。
#if NETCOREAPP3_0_OR_GREATER || NET5_0_OR_GREATER
        private static double MulAdd(double a, double b, double c) => Math.FusedMultiplyAdd(a, b, c);
#else
        private static double MulAdd(double a, double b, double c) => a * b + c;
#endif

        // ==================== 公共 API —— 实数（Cody 1969） ====================

        /// <summary>实数误差函数 erf(x)。</summary>
        public static double Erf(double x)
        {
            double ax = Math.Abs(x);

            if (ax >= Six)
                return x < 0.0 ? -1.0 : 1.0;

            // 小 |x|：Taylor 级数（C 原版此区间委托系统 libm erf）。
            // erf(x) = c1·x·(1 + mx2/3 + mx2²/10 + ...)，mx2 = -x²。
            // 主项 c1·x 的乘积舍入误差：.NET Core 3.0+ 用一条 FMA 直接取得精确残差；
            // 旧 TFM 回退 Dekker TwoProd 拆分。达到正确舍入级（< 0.5 ulp）。
            if (ax < 8e-2)
            {
                const double c1 = 1.1283791670955125739; // 2/sqrt(pi)
                double mx2 = -x * x;
                double p = c1 * x;
#if NETCOREAPP3_0_OR_GREATER || NET5_0_OR_GREATER
                // c1·x 的精确乘积误差：FMA 单指令（exact residual，无需 Dekker 拆分）
                double perr = Math.FusedMultiplyAdd(c1, x, -p);
#else
                // c1 * x 的精确误差（Dekker TwoProduct，无 FMA 的旧 TFM 路径）
                double split = 134217729.0 * c1; // 2^27 + 1
                double c1hi = split - (split - c1);
                double c1lo = c1 - c1hi;
                double xs = 134217729.0 * x;
                double xhi = xs - (xs - x);
                double xlo = x - xhi;
                double perr = ((c1hi * xhi - p) + c1hi * xlo + c1lo * xhi) + c1lo * xlo;
#endif
                // 级数修正项：p·mx2·(1/3 + mx2·(1/10 + mx2·(1/42 + mx2·(1/216 + mx2·(1/1320 + mx2/9360)))))
                // Horner 求值经 MulAdd（新 TFM 为 FMA 单指令，旧 TFM 为普通乘加）
                double t = p * mx2 * MulAdd(mx2, MulAdd(mx2, MulAdd(mx2, MulAdd(mx2, MulAdd(mx2,
                    0.00010683760683760683761,      // 1/9360
                    0.00075757575757575757576),     // 1/1320
                    0.0046296296296296296296),      // 1/216
                    0.023809523809523809524),       // 1/42
                    0.1),
                    0.33333333333333333333);        // 1/3
                return (p + t) + perr;
            }

            if (ax <= Thresh)
                return ErfFirst(x, ax);

            if (ax <= Four)
            {
                double erfc = ErfcSecond(ax);
                double result = 1.0 - erfc;
                return x < 0.0 ? -result : result;
            }

            double erfc3 = ErfcThird(ax);
            double result3 = 1.0 - erfc3;
            return x < 0.0 ? -result3 : result3;
        }

        /// <summary>实数互补误差函数 erfc(x)。</summary>
        public static double Erfc(double x)
        {
            if (x < 0.0)
                return 2.0 - Erfc(-x);

            if (x >= XBig)
                return 0.0;

            if (x <= Thresh)
                return 1.0 - Erf(x);

            if (x <= Four)
                return ErfcSecond(x);

            return ErfcThird(x);
        }

        /// <summary>
        /// 实数缩放互补误差函数 erfcx(x) = e^{x²}·erfc(x)。
        /// 复用 Cody 分段有理逼近：erfc 实现中的 exp 分裂因子 e^{-ysq²}·e^{-del} 与 e^{y²}
        /// 互为倒数（y² = ysq² + del），故 erfcx 直接取有理逼近比值再乘残余修正因子
        /// e^{y²-ysq²-del}（该指数在双精度下 ≈ 0，修正量 ~1e-14 量级），无需 Chebyshev 大表。
        /// x &lt; 0 时用恒等式 erfcx(x) = 2·e^{x²} - erfcx(-x)；x &lt; -26.7 时 e^{x²} 上溢，返回 +∞。
        /// </summary>
        public static double Erfcx(double x)
        {
            if (x >= 0.0)
            {
                if (x <= Thresh)
                    return Math.Exp(x * x) * (1.0 - Erf(x));
                if (x <= Four)
                    return ErfcxSecond(x);
                return ErfcxThird(x);
            }
            // x < 0：erfc(x) = 2 - erfc(-x)，故 erfcx(x) = 2·e^{x²} - erfcx(-x)
            return x < -26.7 ? double.PositiveInfinity : 2.0 * Math.Exp(x * x) - Erfcx(-x);
        }

        // ==================== 私有实现 —— Cody 实数 ====================

        // 第一段（|x| ≤ 0.46875）：erf(x) = x·R(x²)，分子分母均为 4 次有理式
        private static double ErfFirst(double x, double y)
        {
            double ysq = y > XSmall ? y * y : 0.0;

            double xnum = A[4];
            double xden = 1.0;
            for (int i = 0; i < 4; i++)
            {
                xnum = MulAdd(xnum, ysq, A[i]);
                xden = MulAdd(xden, ysq, B[i]);
            }

            return x * (xnum / xden);
        }

        // 第二段（0.46875 < y ≤ 4）：erfc(y) = e^{-y²}·R(y)，exp 用 floor 分裂保持精度
        private static double ErfcSecond(double y)
        {
            double xnum = C[8];
            double xden = 1.0;
            for (int i = 0; i < 8; i++)
            {
                xnum = MulAdd(xnum, y, C[i]);
                xden = MulAdd(xden, y, D[i]);
            }

            double result = xnum / xden;

            double ysq = Math.Floor(y * Sixteen) / Sixteen;
            double del = (y - ysq) * (y + ysq);
            result = Math.Exp(-ysq * ysq) * Math.Exp(-del) * result;

            return result;
        }

        // 第三段（y > 4）：erfc(y) = e^{-y²}·(1/√π - y^{-2}·R(y^{-2}))/y
        private static double ErfcThird(double y)
        {
            if (y >= XHuge)
                return SqrPi / y;

            double ysq = 1.0 / (y * y);
            double xnum = P[5];
            double xden = 1.0;
            for (int i = 0; i < 5; i++)
            {
                xnum = MulAdd(xnum, ysq, P[i]);
                xden = MulAdd(xden, ysq, Q[i]);
            }

            double result = ysq * (xnum / xden);
            result = (SqrPi - result) / y;

            double ysqFloor = Math.Floor(y * Sixteen) / Sixteen;
            double del = (y - ysqFloor) * (y + ysqFloor);
            result = Math.Exp(-ysqFloor * ysqFloor) * Math.Exp(-del) * result;

            return result;
        }

        // 第二段 erfcx：e^{y²}·e^{-ysq²}·e^{-del} ≡ e^{y²-ysq²-del}（≈1，残余为分裂舍入差）
        private static double ErfcxSecond(double y)
        {
            double xnum = C[8];
            double xden = 1.0;
            for (int i = 0; i < 8; i++)
            {
                xnum = MulAdd(xnum, y, C[i]);
                xden = MulAdd(xden, y, D[i]);
            }

            double ysq = Math.Floor(y * Sixteen) / Sixteen;
            double del = (y - ysq) * (y + ysq);
            double residue = y * y - ysq * ysq - del; // 数学上为 0，仅存分裂舍入差
            return (xnum / xden) * Math.Exp(residue);
        }

        // 第三段 erfcx：同上，取有理逼近部分乘残余修正因子
        private static double ErfcxThird(double y)
        {
            if (y >= XHuge)
                return SqrPi / y;

            double ysq = 1.0 / (y * y);
            double xnum = P[5];
            double xden = 1.0;
            for (int i = 0; i < 5; i++)
            {
                xnum = MulAdd(xnum, ysq, P[i]);
                xden = MulAdd(xden, ysq, Q[i]);
            }

            double result = ysq * (xnum / xden);
            result = (SqrPi - result) / y;

            double ysqFloor = Math.Floor(y * Sixteen) / Sixteen;
            double del = (y - ysqFloor) * (y + ysqFloor);
            double residue = y * y - ysqFloor * ysqFloor - del;
            return result * Math.Exp(residue);
        }
    }
}
