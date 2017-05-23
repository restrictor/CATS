nb = c(73,78,72,72,72,81)
svm = c(86,82,86,87,86,84)
nn = c(88,84,89,89,89,86)
mr = c(81,84,86,86,86,87)
rf = c(86,84,85,85,81,87)

SIMPLE = data.frame(nb,svm,nn,mr,rf)
SIMPLE_mean = apply(SIMPLE, 2, mean)
SIMPLE_sd = apply(SIMPLE, 2, sd)
SIMPLE_mean
SIMPLE_sd

RFSVM = c(86,88,90,86,88,85,89)
RF2SVM = c(93,90,93,93,93,92,93)
NNET2SVM = c(91,91,93,92,92,91,91)
SVM2SVM = c(91,90,93,92,92,91,91)
SVM2SVMRF = c(91,87,91,91,91,90,88)

EASY = data.frame(RFSVM,RF2SVM,NNET2SVM,SVM2SVM,SVM2SVMRF)
EASY_mean = apply(EASY, 2, mean)
EASY_sd = apply(EASY, 2, sd)
EASY_mean
EASY_sd
