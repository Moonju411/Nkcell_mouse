1. Anaconda설치
bash /home/shared/program/anaconda/linux/Anaconda3-2020.11-Linux-x86_64.sh

2. Anaconda environment 설정
conda create -n monocle3 r-essentials r-base #monocle3라는 이름의 conda envronmnet만들고, r 설치
conda activate monocle3 #monocle3라는 이름의 conda envronmnet를 로딩
conda install -c bioconda r-monocle3 #monocle3라는 이름의 conda envronmnet에 monocle3 설치
conda install -c conda-forge r-seurat #monocle3라는 이름의 conda envronmnet에 seurat 설치 ( r 버전때문에 seurat 버전이 달라짐)

3. 주피터 노트북 설치/설정
jupyter kernelspec list으로 ir 뜨는 지 확인
![화면 캡처 2021-04-07 195820](https://user-images.githubusercontent.com/42495757/113856052-b24db080-97db-11eb-86ed-98364819beee.png)

4. 주피터 노트북에서
install.packages("remotes")
remotes::install_github('satijalab/seurat-wrappers’)

