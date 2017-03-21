///////////////////////////////////////////////////////////////////////////////
//
// Ellipsoid fitting for magnetic sensor hard-iron and soft-iron calibration
// Based on Q. Li. algorithm detailed here:
// https://sites.google.com/site/sailboatinstruments1/step-1
//
// For compiling you will need Math.NET package from NuGet
// TabSize: 2
///////////////////////////////////////////////////////////////////////////////
using MathNet.Numerics.Data.Text;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using System;
using System.Numerics;

namespace MagCal
{
	class Program
	{
		static void Main(string[] args)
		{
			string filename = @"mag.txt";
			int i;

			// load measurement data
			Matrix<double> data = DelimitedReader.Read<double>(filename, false, "\t", true);

			// create D 10xnumdata
			Matrix<double> D = Matrix<double>.Build.Dense(10, data.RowCount);

			for (i = 0; i < data.RowCount; i++)
			{
				D[0, i] = data[i, 0] * data[i, 0];
				D[1, i] = data[i, 1] * data[i, 1];
				D[2, i] = data[i, 2] * data[i, 2];
				D[3, i] = 2.0 * data[i, 1] * data[i, 2];
				D[4, i] = 2.0 * data[i, 0] * data[i, 2];
				D[5, i] = 2.0 * data[i, 0] * data[i, 1];
				D[6, i] = 2.0 * data[i, 0];
				D[7, i] = 2.0 * data[i, 1];
				D[8, i] = 2.0 * data[i, 2];
				D[9, i] = 1.0;
			}

			// create S 10x10
			Matrix<double> S = D * D.Transpose();

			// create sub matrices of S
			//
			//      -                   -
			// S = | S11(6x6)   S12(6x4) |
			//     | S12_t(4x6) S22(4x4) |
			//      -                   -
			//

			// create S11   6x6
			Matrix<double> S11 = Matrix<double>.Build.Dense(6, 6);
			S11.SetSubMatrix(0, 0, 6, 0, 0, 6, S);


			// create S12  6x4
			Matrix<double> S12 = Matrix<double>.Build.Dense(6, 4);
			S12.SetSubMatrix(0, 0, 6, 0, 6, 4, S);

			// create S12t  4x6
			Matrix<double> S12t = Matrix<double>.Build.Dense(4, 6);
			S12t.SetSubMatrix(0, 6, 4, 0, 0, 6, S);

			// create S22 4x4
			Matrix<double> S22 = Matrix<double>.Build.Dense(4, 4);
			S22.SetSubMatrix(0, 6, 4, 0, 6, 4, S);

			// calculate pseudo inverse of S22
			Matrix<double> S22_1 = S22.PseudoInverse();

			// calculate SS = S11 - S12 * S22_1 * S12t
			Matrix<double> SS = S11 - S12 * S22_1 * S12t;

			// Create constraint matrix C
			Matrix<double> Co = Matrix<double>.Build.DenseOfArray(new double[,] {
			{ -1.0,  1.0,  1.0,  0.0,  0.0,  0.0 },
			{  1.0, -1.0,  1.0,  0.0,  0.0,  0.0 },
			{  1.0,  1.0, -1.0,  0.0,  0.0,  0.0 },
			{  0.0,  0.0,  0.0, -4.0,  0.0,  0.0 },
			{  0.0,  0.0,  0.0,  0.0, -4.0,  0.0 },
			{  0.0,  0.0,  0.0,  0.0,  0.0, -4.0 } });

			Matrix<double> C = Co.Inverse();

			// Calculate E = C * SS
			Matrix<double> E = C * SS;

			// calculate eigenvalues wr(6x1) and eigenvectors vr(6x6) of matrix E
			Evd<double> eigen = E.Evd();
			Vector<Complex> wr = eigen.EigenValues;
			Matrix<double> vr = eigen.EigenVectors;

			// find the zero based position of the only positive eigenvalue. The associated eigenvector will be in the corresponding column of matrix vr 
			int index = 0;
			double maxval = wr[0].Real;
			for (i = 1; i < 6; i++)
			{
				if (wr[i].Real > maxval)
				{
					maxval = wr[i].Real;
					index = i;
				}
			}

			// Extract the associated eigenvector v1
			Vector<double> v1 = vr.Column(index);

			// check sign of eigenvector v1
			if (v1[0] < 0.0)
			{
				v1[0] = -v1[0];
				v1[1] = -v1[1];
				v1[2] = -v1[2];
				v1[3] = -v1[3];
				v1[4] = -v1[4];
				v1[5] = -v1[5];
			}

			// Calculate v2 = S22a * v1 
			Vector<double> v2 = (S22_1 * S12t) * v1;

			// calculate v
			Vector<double> v = Vector<double>.Build.Dense(10);

			v[0] = v1[0];
			v[1] = v1[1];
			v[2] = v1[2];
			v[3] = v1[3];
			v[4] = v1[4];
			v[5] = v1[5];
			v[6] = -v2[0];
			v[7] = -v2[1];
			v[8] = -v2[2];
			v[9] = -v2[3];

			//At this point, we have found the general equation of the fitted ellipsoid:
			// Ax² + By² + Cz² + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hg + 2Iz + J = 0

			// where:
			// A = v[0] - term in x2
			// B = v[1] - term in y2
			// C = v[2] - term in z2
			// D = v[5] - term in xy
			// E = v[4] - term in xz
			// F = v[3] - term in yz
			// G = v[6] - term in x
			// H = v[7] - term in y
			// I = v[8] - term in z
			// J = v[9] - constant term

			// If we define
			//      -     -          - -
			//     | A D E |        | G |
			// Q = | D B F |    U = | H |
			//     | E F C |        | I |
			//      -     -          - -
			//                                                                        
			// then the center of the ellipsoid can be calculated as the vector B = -Qˉ¹ * U.
			// The center of the ellipsoid represents the combined bias.

			Matrix<double> Q = Matrix<double>.Build.Dense(3,3);

			Q[0,0] = v[0]; // A
			Q[0,1] = v[5]; // D
			Q[0,2] = v[4]; // E
			Q[1,0] = v[5]; // D
			Q[1,1] = v[1]; // B
			Q[1,2] = v[3]; // F
			Q[2,0] = v[4]; // E
			Q[2,1] = v[3]; // F
			Q[2,2] = v[2]; // C

			Vector<double> U = Vector<double>.Build.Dense(3);

			U[0] = v[6]; // G
			U[1] = v[7]; // H
			U[2] = v[8]; // I

			// Calculate matrix Q_1, the inverse of matrix Q
			Matrix<double> Q_1 = Q.Inverse();

			// Calculate B = Q_1 * U   ( 3x1 = 3x3 * 3x1)
			Vector<double> B = Q_1 * U;

			// Calculate combined bias
			B[0] = -B[0];     // x-axis combined bias
			B[1] = -B[1];     // y-axis combined bias
			B[2] = -B[2];     // z-axis combined bias

			Console.WriteLine("Combined bias:");
			Console.WriteLine(B.ToString());

			//            -1
			// Calculate A
			//
			//  -1         Hm                 1/2
			// A    = -------------------- * Q
			//        sqrt(Bt * Q * B - J)


			// Calculate btqb = BT * Q * B
			double btqb = B * Q * B;

			// Calculate hmb = sqrt(btqb - J).
			double J = v[9];
			double hmb = Math.Sqrt(btqb - J);

			// Calculate SQ, the square root of matrix Q
			eigen = Q.Evd();
			wr = eigen.EigenValues;
			vr = eigen.EigenVectors;

			// normalize eigenvectors
			double norm1 = Math.Sqrt(vr[0,0] * vr[0,0] + vr[0,1] * vr[0,1] + vr[0,2] * vr[0,2]);
			vr[0,0] /= norm1;
			vr[0,1] /= norm1;
			vr[0,2] /= norm1;
			double norm2 = Math.Sqrt(vr[1,0] * vr[1,0] + vr[1,1] * vr[1,1] + vr[1,2] * vr[1,2]);
			vr[1,0] /= norm2;
			vr[1,1] /= norm2;
			vr[1,2] /= norm2;
			double norm3 = Math.Sqrt(vr[2,0] * vr[2,0] + vr[2,1] * vr[2,1] + vr[2,2] * vr[2,2]);
			vr[2,0] /= norm3;
			vr[2,1] /= norm3;
			vr[2,2] /= norm3;

			Matrix<double> Dz = Matrix<double>.Build.Dense(3, 3);
			Dz[0, 0] = Math.Sqrt(wr[0].Real);
			Dz[1, 1] = Math.Sqrt(wr[1].Real);
			Dz[2, 2] = Math.Sqrt(wr[2].Real);

			Matrix<double> SQ = (vr * Dz) * vr.Transpose();

			double hm = 0.569;

			Matrix<double> A_1 = SQ * hm / hmb;

			Console.WriteLine("A-1 matrix:");
			Console.WriteLine(A_1.ToString());

			// Calculate A to permit comparison with MagCal 
			Matrix<double> A = A_1.Inverse();

			Console.WriteLine("A matrix:");
			Console.WriteLine(A.ToString());
		}
	}
}
