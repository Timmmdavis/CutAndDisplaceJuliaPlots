function CptReader(Fid)
#Lazy way of loading CPT files from: http://soliton.vm.bytemark.co.uk/pub/cpt-city/km/tn/cool-warm-d12.png.index.html


#Fid=raw"C:\Users\timmm\Documents\GitHub\CutAndDisplace\SharedFunctions\Ccool-warm2.cpt"
#(cmap2)=CptReader(Fid)

#Open the file
#using DelimitedFiles
TextFile=readdlm(Fid);

#Size of the array
n=length(TextFile[:,1]);

#Find row index's of values
FirstRow=0; 
LastRow=0;
for i=1:n
	if TextFile[i,1]!="#"
		FirstRow=i;
		break
	end
end	
for i=1:n
	if TextFile[end-i,1]=="B"
		LastRow=n-i-1;
		break
	end
end	
#Extract from data (if you get indexing errors check the last 3 rows are something like:
#B	0	0	0
#F	255	255	255
#N	255	0	0 )
#and the first few:
#	cpt file created by: makecpt -T-6.0/6.0/1 -Ccool-warm.cpt
#COLOR_MODEL = RGB
#
CPT=TextFile[FirstRow:LastRow,2:4]
#CPT=[CPT ones(length(CPT[:,1]))] #Add alpha channel


return CPT
end