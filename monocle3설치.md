1. Anaconda설치  
`bash /home/shared/program/anaconda/linux/Anaconda3-2020.11-Linux-x86_64.sh`

2. Anaconda environment 설정  
`conda create -n monocle3 r-essentials r-base` #monocle3라는 이름의 conda envronmnet만들고, r 설치  
`conda activate monocle3` #monocle3라는 이름의 conda envronmnet를 로딩  
`conda install -c bioconda r-monocle3` #monocle3라는 이름의 conda envronmnet에 monocle3 설치  
`conda install -c conda-forge r-seurat` #monocle3라는 이름의 conda envronmnet에 seurat 설치 ( r 버전때문에 seurat 버전이 달라짐)  

2-1. 미리 설정된 .yaml파일이 있는 경우
`conda env create -n {env name} --file {environment.yaml}`

3. 주피터 노트북 설치/설정
`jupyter kernelspec list`으로 ir 뜨는 지 확인. 
![화면 캡처 2021-04-07 195820](https://user-images.githubusercontent.com/42495757/113856052-b24db080-97db-11eb-86ed-98364819beee.png)

4. 주피터 노트북 접속  
`jupyter notebook --ip=172.27.30.99`
빨간줄 복사에서 웹브라우저로 접속하기
![image](https://user-images.githubusercontent.com/42495757/113856474-41f35f00-97dc-11eb-827f-8703c331094a.png)

5. 주피터 노트북에서
`install.packages("remotes")`
`remotes::install_github('satijalab/seurat-wrappers’)` #seurat이랑 monocle3 같이쓰기

