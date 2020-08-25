package net.dedekind.blas;

import java.util.Arrays;

public class ComplexGemmTest {

    public ComplexGemmTest() {
        // expected result
        // (0, 0): 6.575 3.45
        // (0, 1): 8.5375 -0.1499999999999999
        // (0, 2): 8.375 4.75
        // (1, 0): 28.225 10.35
        // (1, 1): 34.9875 4.550000000000001
        // (1, 2): 36.375 20.25
    }

    public static void main(String[] args) {
        // zgemm
        // 2 x 3 matrix C
        double[] c = new double[12];

        // 2 x 2 matrix A
        // 0.0 + 0.0i, 1.0 + 0.25i
        // 2.0 + 0.5i, 3.0 + 0.75i
        double[] a = new double[8];
        //@formatter:off
        a[0] = 0.0; a[1] = 0.0;
        a[2] = 2.0; a[3] = 0.5;
        a[4] = 1.0; a[5] = 0.25;
        a[6] = 3.0; a[7] = 0.75;
        //@formatter:on

        // 2 x 3 matrix B
        // 4.0 - 1.0i, 5.0 + 1.25i, 6.0 + 1.5i
        // 7.0 + 1.7i, 8.0 - 2.15i, 9.0 + 2.5i
        double[] b = new double[12];
        //@formatter:off
        b[0] = 4.0; b[1] = -1.0;
        b[2] = 7.0; b[3] = 1.7;
        b[4] = 5.0; b[5] = 1.25;
        b[6] = 8.0; b[7] = -2.15;
        b[8] = 6.0; b[9] = 1.5;
        b[10] = 9.0; b[11] = 2.5;
        //@formatter:on

        int numRowsC = 2;
        int numColsC = 3;

        int numRowsA = 2;
        int numColsA = 2;

        int numRowsB = 2;
        int numColsB = 3;

        double alphar = 1.0;
        double alphai = 0.0;
        double betar = 1.0;
        double betai = 0.0;

        System.out.println("A: " + Arrays.toString(a));
        Blas blas = Blas.getInstance();
        blas.zgemm(Trans.N, Trans.N, numRowsC, numColsC, numColsA, alphar, alphai, a, Math.max(1, numRowsA), b,
                Math.max(1, numRowsB), betar, betai, c, Math.max(1, numRowsC));

        System.out.println("C: " + Arrays.toString(c));
        System.out.println("A: " + Arrays.toString(a));

        // cgemm
        float alpharf = (float) alphar;
        float alphaif = (float) alphai;
        float betarf = (float) betar;
        float betaif = (float) betai;
        float[] c2 = new float[12];
        float[] a2 = new float[8];
        float[] b2 = new float[12];
        for (int i = 0; i < a2.length; ++i) {
            a2[i] = (float) a[i];
        }
        for (int i = 0; i < b2.length; ++i) {
            b2[i] = (float) b[i];
        }

        System.out.println(" --- ");
        System.out.println("A: " + Arrays.toString(a2));
        blas.cgemm(Trans.N, Trans.N, numRowsC, numColsC, numColsA, alpharf, alphaif, a2, Math.max(1, numRowsA), b2,
                Math.max(1, numRowsB), betarf, betaif, c2, Math.max(1, numRowsC));

        System.out.println("C: " + Arrays.toString(c2));
        System.out.println("A: " + Arrays.toString(a2));
    }
}
