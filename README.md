# CutAndDisplaceJuliaPlots


Example:
#Run: CutAndDisplaceJulia: test_TDvsEshelbyPennyCrack.
using CutAndDisplaceJuliaPlots
Fid=CutAndDisplaceJulia.LoadData(CutAndDisplaceJuliaPlots,"Ccool-warm2.cpt")
(cmap2)=CutAndDisplaceJuliaPlots.CptReader(Fid);
cmap2=cmap2./255; #RGB2Colour
using CutAndDisplaceJuliaPlots
CutAndDisplaceJuliaPlots.DrawMeshMakie(TotalShearing,"TotalShearing",P1,P2,P3,cmap2)
test=1 #So it all pastes in
