# R
##清除环境中所有对象
rm(list=ls())

##安装包
library(pracma)
library(drc)
library("TH.data")

##设置药物浓度
conc <- c(0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 200)

##Vehicle response（溶剂对照组响应值）
Control <- c(
  86.42959585,
  82.67949039,
  62.39348686,
  70.86813462,
  67.23614272,
  69.53935710,
  48.83995611,
  25.66017041,
  0.41339745,
  1.06302202
)

##Modulator response（调节剂组响应值）
modulator <- c(
  99.54413673,
  97.05760464,
  93.03771239,
  102.36220474,
  92.54040615,
  87.89888106,
  71.73642768,
  30.50145048,
  0.49730626,
  1.36759221
)

##浓度-响应值数据整合
vehicle_data <- data.frame(
  conc=conc,
  response=Control
)

modulator_data <- data.frame(
  conc=conc,
  response=modulator
)

##4PL拟合（量效关系）
fit_vehicle <- drm(
  response ~ conc,
  data=vehicle_data,
  fct=LL.4()
)

fit_modulator <- drm(
  response ~ conc,
  data=modulator_data,
  fct=LL.4()
)

##生成拟合曲线
new_conc <- 10^seq(
  log10(min(conc)),
  log10(max(conc)),
  length.out=500
)

pred_vehicle <- predict(
  fit_vehicle,
  newdata=data.frame(conc=new_conc)
)

pred_modulator <- predict(
  fit_modulator,
  newdata=data.frame(conc=new_conc)
)

##计算拟合曲线AUC
AUC_vehicle <- trapz(
  log10(new_conc),
  pred_vehicle
)

AUC_modulator <- trapz(
  log10(new_conc),
  pred_modulator
)

##DAUCn计算
DAUCn <- (AUC_modulator-AUC_vehicle)/
  AUC_vehicle

print(DAUCn)

##绘制溶剂对照组曲线
plot(
  new_conc,
  pred_vehicle,
  type="l",
  log="x",
  col="#5A97D0",
  lwd=3,
  ylim=c(0,110),
  xlab="Concentration",
  ylab="Response (%)",
  main=paste0(
    "Dose-response curve\nDAUCn = ",
    round(DAUCn,3)
  ),
  font.lab = 2,
  font.main = 2,
  font.axis = 2,
  las = 1,
)

##添加溶剂对照组实验点
points(
  conc,
  Control,
  col="#5A97D0",
  pch=16
)

##绘制调节剂组曲线
lines(
  new_conc,
  pred_modulator,
  col="#F0A0A0",
  lwd=3
)

##添加调节剂组实验点
points(
  conc,
  modulator,
  col="#F0A0A0",
  pch=16
)

##添加图注
legend(
  "topright",
  legend=c(
    "sgAAVS1",
    "sgBRCA1#1"
  ),
  col=c("#5A97D0","#F0A0A0"),
  lwd=3,
  pch=16,
  text.font = 2,
  cex = 0.8
)

##颜色填充
polygon(
  c(new_conc, rev(new_conc)),
  c(pred_vehicle, rev(pred_modulator)),
  col=rgb(1,0,0,0.25),
  border=NA
)
