package model;

public class FastAtan {

    private static final int ATAN_2_BITS = 7;

    private static final int ATAN_2_BITS2 = ATAN_2_BITS << 1;
    private static final int ATAN_2_MASK = ~(-1 << ATAN_2_BITS2);
    private static final int ATAN_2_COUNT = ATAN_2_MASK + 1;
    private static final int ATAN_2_DIM = (int) Math.sqrt(ATAN_2_COUNT);

    private static final double INV_ATAN_2_DIM_MINUS_1 = 1.0f / (ATAN_2_DIM - 1);

    private static final double[] atan2 = new double[ATAN_2_COUNT];



    static {
        for (int i = 0; i < ATAN_2_DIM; i++) {
            for (int j = 0; j < ATAN_2_DIM; j++) {
                double x0 = (double) i / ATAN_2_DIM;
                double y0 = (double) j / ATAN_2_DIM;

                atan2[j * ATAN_2_DIM + i] = Math.atan2(y0, x0);
            }
        }
    }


    /**
     * ATAN_2
     */

    public static double atan2(double y, double x) {
        double add, mul;

        if (x < 0.0f) {
            if (y < 0.0f) {
                x = -x;
                y = -y;

                mul = 1.0f;
            } else {
                x = -x;
                mul = -1.0f;
            }

            add = -3.141592653f;
        } else {
            if (y < 0.0f) {
                y = -y;
                mul = -1.0f;
            } else {
                mul = 1.0f;
            }

            add = 0.0f;
        }

        double invDiv = 1.0f / ((Math.max(x, y)) * INV_ATAN_2_DIM_MINUS_1);

        int xi = (int) (x * invDiv);
        int yi = (int) (y * invDiv);

        return (atan2[yi * ATAN_2_DIM + xi] + add) * mul;
    }
}
