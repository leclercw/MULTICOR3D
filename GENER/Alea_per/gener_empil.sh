disx=100 # discretisation
disy=100 # discretisation
disz=100 # discretisation
dimx=100 # dimension de la boite
dimy=100 # dimension de la boite
dimz=100 # dimension de la boite
nbp=6000000 # nombre de particules
fs=0.64 # fraction surfacique de particules
nmap=1 # nmap
sep="_"

export OMP_STACKSIZE=10g

for ((k=0 ; $nmap - $k; k++))
do 
./gener3d_coh $disx $disy $disz $dimx $dimy $dimz $nbp $fs
filenameo="test${nbp}_${k}.epl"
mv pristimap $filenameo
 
done


