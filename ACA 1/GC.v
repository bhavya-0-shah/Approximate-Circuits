module GC #(parameter valency=4) (
	output GG, 
	input [valency-1 : 0] g, 
	input [valency-1 : 1] p
);
	wire [valency-1 : 0] wr, gg;

	assign gg[0]=g[0];
	
	genvar k;

	generate
		for(k=0; k<valency-1; k=k+1) begin
			and (wr[k], p[k+1], gg[k]);
			or (gg[k+1], wr[k], g[k+1]);                   // For valency 4, G[3:0]=G3+P3G2+P3P2G1+P3P2P1G0     
		end
	endgenerate	

	assign GG=gg[valency-1];
endmodule