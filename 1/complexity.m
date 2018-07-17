
Nsat = 37;
Nf = 41*(1:16);
M = 1e4;
NT = 1e3*(1:16);
figure();
[DFC_multi,DFC_add] = arrayfun(@(x,y) DFC_comp(x,Nsat,y,M),NT,Nf);
plot(NT,DFC_multi);
[NCH_multi,NCH_add] = arrayfun(@(x) NCH_comp(x,Nsat,41,M),10);
hold on
plot(NT,repmat(NCH_multi,[1,16]));
figure();
plot(NT,DFC_add);
hold on
plot(NT,repmat(NCH_add,[1,16]));
function [multiNum,addNum ] = DFC_comp(NT,Nsat,Nf,M)
multiNum = M*(NT+4*Nf*Nsat+4*Nf*Nsat*log2(M));
addNum = M*Nf*(NT+2*Nsat+8*Nsat*log2(M));
end
function [multiNum,addNum ] = NCH_comp(Nc,Nsat,Nf,M)
multiNum = M*Nf*(2*Nc+4*Nsat*Nc+4*Nsat*Nc*log2(M));
addNum = M*Nf*(4*Nsat*Nc+8*Nc*Nsat*log2(M));
end