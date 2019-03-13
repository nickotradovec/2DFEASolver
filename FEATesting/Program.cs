using System;
using System.Collections;
using System.Collections.Generic;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace FEATesting
{
    class Program
    {
        const double E = 70000000000; // Young's Modulus (Pa)
        const double nu = 0.33; // Poisson's ratio (unitless)

        static void Main(string[] args)
        {
            // Proof of concept. See example 8.1 at:
            // http://www2.eng.cam.ac.uk/~jl305/3D7/slid.pdf
            // We will gradually increase generalization of this problem.
            
            Matrix<double> DPlaneStress = DenseMatrix.OfArray(new double[,] {
                {1, nu, 0},
                {nu, 1, 0},
                {0, 0, (1-nu)/2}}) * (E / (1 - Math.Pow(nu, 2)));

            double[,] nodes = {{0, 0}, {2, 0}, {2, 1}, {0, 1}}; // (x1, y1), (x2, y2)...
            
            int[] element1 = {0, 2, 3};
            int[] element2 = {0, 1, 2};        
            int[][] elements = {element1, element2}; // Element1[node1, node2, node3], ...

            // Compute Strain-Displacement Matrices
            Matrix<double>[] B = new Matrix<double>[elements.Length];
            for(int i=0; i<elements.Length; i++) {
                B[i] = BElement(nodes, elements[i], 1D); // Area should be computed
            }

            // Compute Local Stiffnesses
            Matrix<double>[] K = new Matrix<double>[elements.Length];
            for(int i=0; i<elements.Length; i++) {
                K[i] = (B[i].Transpose()) * DPlaneStress * B[i] * 1D;
            }

            // Assemble Global Stiffness Matrix
            // 2D, so each node as 2 DOF
            // DOF nodes will be represented by the node index * 2 (+1 for y)
            // TODO: Figure out a better way to assemble this matrix
            Matrix<double> KGlobal = new DenseMatrix(nodes.GetLength(0) * 2, nodes.GetLength(0) * 2);
            for(int i=0; i<K.Length; i++) {

                int[] mapping = Mapping(elements[i]);

                for(int j=0; j<K[i].RowCount; j++) {
                    for(int k=0; k<K[i].ColumnCount; k++) {
                        
                        KGlobal[mapping[j], mapping[k]] += K[i][j, k];
                    }
                }
            }

            // Eliminate rows and columns as per boundary conditions   
            // TODO: Here we kind of lose track of what indeces represent what. This should be improved.  
            int[] noDoF = {7, 6, 1, 0}; // Specify in reverse order for clarity on indexes
            foreach (var val in noDoF) { 
                KGlobal = KGlobal.RemoveColumn(val); 
                KGlobal = KGlobal.RemoveRow(val); 
            }

            // Apply External Force and Solve
            var externalForce = new DenseMatrix(4, 1);
            externalForce[0, 0] = 5;
            externalForce[1, 0] = 10;
            externalForce[2, 0] = 5;
            externalForce[3, 0] = 10;

            Matrix<double> KDisplacement = KGlobal.Solve(externalForce);
        }

        static Matrix<double> BElement(double[,] nodes, int[] element, double A)
        {
            var rtn = new DenseMatrix(3, 6); // rows, columns
            rtn[0, 0] = (nodes[element[1], 1] - nodes[element[2], 1]);
            rtn[0, 2] = (nodes[element[2], 1] - nodes[element[0], 1]);
            rtn[0, 4] = (nodes[element[0], 1] - nodes[element[1], 1]);

            rtn[1, 1] = (nodes[element[2], 0] - nodes[element[1], 0]);
            rtn[1, 3] = (nodes[element[0], 0] - nodes[element[2], 0]);
            rtn[1, 5] = (nodes[element[1], 0] - nodes[element[0], 0]);

            rtn[2, 0] = (nodes[element[2], 0] - nodes[element[1], 0]);
            rtn[2, 1] = (nodes[element[1], 1] - nodes[element[2], 1]);
            rtn[2, 2] = (nodes[element[0], 0] - nodes[element[2], 0]);
            rtn[2, 3] = (nodes[element[2], 1] - nodes[element[0], 1]);
            rtn[2, 4] = (nodes[element[1], 0] - nodes[element[0], 0]);
            rtn[2, 5] = (nodes[element[0], 1] - nodes[element[1], 1]);

            return rtn * (1 / (2 * A));
        }

        static int[] Mapping(int[] element) {

            var rtn = new int[element.Length * 2];
            for(int i=0; i<element.Length; i++) {
                rtn[2 * i] = 2 * element[i];
                rtn[(2 * i) + 1] = (2 * element[i]) + 1;
            }
            return rtn;
        }
    }
}