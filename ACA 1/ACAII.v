`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer:        Pravin Mane
// 
// Create Date:    
// Design Name: 
// Module Name:   ACA_II
// Project Name:  Accuracy Configurabe Adder - II (ACA - II)
// Target Devices: 
// Tool versions: 
// Description: 
//
// Dependencies: 
//
// Revision: 
// Revision 0.01 - File Created
// Additional Comments: 
//
//////////////////////////////////////////////////////////////////////////////////
module ACAII #(parameter N = 16, R = 4, P = 4)(
                output [N:1] SUM, 
                output COUT,
                input [N:1] A, B,
                input CIN);

        wire [(N/R)-1:1] couti;                 //to obtain the carry out produced by each of the sub-adder units
        wire [P:1] notused;                     //prediction bit declaration

                                                /*to calculate the sum of the last sub-adder with 8 bits using CLA Adder.
                                                The last sub-adder produces an (R+P) bit result but the consequent sub adders
                                                produce an R bit result andf the remaining P bits in the sub-adder are the 
                                                prediction bits which predict the carry in to the R bits of the sub-adder which 
                                                determine the result.*/
            CLA_p_v #(R+P,4) CLAINSTL(                               
            SUM[R+P:1],
            couti[1],
            A[R+P:1],
            B[R+P:1],
            CIN
        );

        genvar i; 
        generate
        for (i = 0; i<=(N/R - 2); i = i + 1)  begin:ACA_parameterized   /*to claculate the sum (R bit sum and P prediction bits) produced by the subsequent  
                                                sub-adders.*/
            CLA_p_v #(R+P,4) CLAINST(
                {SUM[N - (i*R) - 1:N - ((i+1)*R)], notused[P:1]},
                    couti[i+2],
                    A[N - (i*R) - 1:N - ((i+2)*R)],
                    B[N - (i*R) - 1:N - ((i+2)*R)],
                    1'b0
                );
            end
        endgenerate
        assign COUT = couti[2];
    endmodule 
