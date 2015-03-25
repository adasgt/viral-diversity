# Shell Script to download all the virus genomes from the genbank 
if [ -d "data" ]; then
	rm data/*
else
	mkdir data
fi
cd data
mkdir tmp
cd tmp
wget -q "ftp://ftp.ncbi.nih.gov/genomes/Viruses/all.fna.tar.gz"
gunzip all.fna.tar.gz
tar -xf all.fna.tar
rm all.fna.tar
for i in `ls`;
do
	cp $i/* ../
done

cd ..
rm -r tmp
ls -l | grep '^-' | awk '{print $10}' > viruses.txt

