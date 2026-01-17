## WT

cd transcript_counts

# WT
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM9438865&format=file&file=GSM9438865%5Fwt%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM9438865&format=file&file=GSM9438865%5Fwt%5F2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"

# DKO
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM9438866&format=file&file=GSM9438866%5Fdko%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM9438866&format=file&file=GSM9438866%5Fdko%5F2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"

# TKO
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM9438867&format=file&file=GSM9438867%5Ftko%5F1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM9438867&format=file&file=GSM9438867%5Ftko%5F2%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5"

# Neural
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3814885&format=file&file=GSM3814885%5Fday0%5Ffiltered%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5"
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3814886&format=file&file=GSM3814886%5Fday4%5Ffiltered%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5"
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3814887&format=file&file=GSM3814887%5Fday6%5Ffiltered%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5"

# ECS
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3814885&format=file&file=GSM3814885%5Fday0%5Ffiltered%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5"
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3814886&format=file&file=GSM3814886%5Fday4%5Ffiltered%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5"
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3814887&format=file&file=GSM3814887%5Fday6%5Ffiltered%5Fgene%5Fbc%5Fmatrices%5Fh5%2Eh5"

# EBs

mkdir EB
cd EB

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE178274&format=file&file=GSE178274%5Fbarcodes%2Etsv%2Egz"
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE178274&format=file&file=GSE178274%5Fcell%5Fcounts%2Emtx%2Egz"
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE178274&format=file&file=GSE178274%5Fcell%5Fmetadata%2Etxt%2Egz"
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE178274&format=file&file=GSE178274%5Fgenes%2Etsv%2Egz"

mv GSE178274_barcodes.tsv.gz barcodes.tsv.gz
mv GSE178274_cell_counts.mtx.gz cell_counts.mtx.gz
mv GSE178274_cell_metadata.txt.gz cell_metadata.txt.gz
mv GSE178274_genes.tsv.gz genes.tsv.gz

cd ../../methylation

mkdir data

# HUES64 DKO
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3662265&format=file&file=GSM3662265%5FHUES64%5FDKO%5FP28%5FWGBS%2Ebed%2Egz"

# HUES64 WT
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4458668&format=file&file=GSM4458668%5FWGBS%5FHUES64%5FWT%5FESCs%2Ebed%2Egz"

# HUES8 TKO
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3618720&format=file&file=GSM3618720%5FHUES8%5FTKO%5FWGBS%2Ebed%2Egz"

# HUES8 WT
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3618718&format=file&file=GSM3618718%5FHUES8%5FWT%5FWGBS%2Ebed%2Egz"