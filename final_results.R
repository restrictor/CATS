nb = c(78,72,72)
svm = c(82,86,87)
nn = c(84,89,89)
mr = c(84,86,86)
rf = c(84,85,85)

SIMPLE = data.frame(nb,svm,nn,mr,rf)
SIMPLE_mean = apply(SIMPLE, 2, mean)
SIMPLE_sd = apply(SIMPLE, 2, sd)
SIMPLE_mean
SIMPLE_sd

RFSVM = c(86,88,90,86,88)
RF2SVM = c(93,90,93,93,93)
NNET2SVM = c(91,91,93,92,92)
SVM2SVM = c(91,90,93,92,92)
SVM2SVMRF = c(91,87,91,91,91)

EASY = data.frame(RFSVM,RF2SVM,NNET2SVM,SVM2SVM,SVM2SVMRF)
EASY_mean = apply(EASY, 2, mean)
EASY_sd = apply(EASY, 2, sd)
EASY_mean
EASY_sd
