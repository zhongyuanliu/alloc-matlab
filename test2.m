x=randn(11,1);
y=randn(size(Aeq,1),1);
z=randn(size(A,1),1);
s=z;
G=H;
g=f;

[ x_stop, y_stop, z_stop, s_stop, k ] = pcQP_gen (x, y, z, s,G, g ,A', b ,Aeq', beq)