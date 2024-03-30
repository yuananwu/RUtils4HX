
# 获取km曲线的文字结果
getKmWordRst = function(kmFit, coxMd, specTmList, kmFo, df){
  kmRst = NULL
  # 输出km曲线中位生存时间
  kmRst[['分组中位生存时间']] = kmFit
  #输出不同时间点的生存概率
  options(max.print = 10000)
  survDetail = survminer::surv_summary(kmFit)
  survDetailRst = as.data.frame(cbind(survDetail$time,as.character(survDetail$strata) ,survDetail$surv,survDetail$lower,survDetail$upper))
  names(survDetailRst) = c('tm','strata','prob','95CILow','95CIUpper')
  getRow = function(df, strata, cutTmPoint){
    rst = df[which(df[, 'strata'] == strata & as.numeric(df[, 'tm']) <= cutTmPoint),]
    if (dim(rst)[1] == 0){ #没有比cutPoint更早的数据
      rst = df[which(df[, 'strata'] == strata),]
      return(rst[1, ])
    }else{# 取得到数, 就取最后一条, 离cutpoint最近的那条
      return(rst[dim(rst)[1], ])
    }
  }
  mergeRst = NULL
  for (fac in unique(survDetailRst$strata)){
    for(tms in specTmList){
      row = getRow(survDetailRst, fac,  tms)
      row[, 'specTm'] = tms
      mergeRst = rbind(mergeRst, row)
    }
  }
  mergeRst = data.frame(mergeRst, row.names = NULL)
  colnames(mergeRst) = c('实际生存时间', '分层标识', '生存概率', '生存概率_95ciLow', '生存概率_95ciUpper', '指定时间点')
  kmRst[['指定时间点的生存概率']] = mergeRst
  kmRst[['km曲线总体对比']] = survival::survdiff(kmFo, df)$pvalue
  kmRst[['km曲线两两对比']] = survminer::pairwise_survdiff(kmFo, df, p.adjust.method = "BH")
  kmRst[['单因素cox得到的HR等指标']] = summary(coxMd)$coefficients
  kmRst[['c-index']] = getCindex95CI(coxMd, df, roundNum = 3)
  return(kmRst)
}

#将txt格式的kmRst整理成csv格式
kmRstReshape = function(kmRst){
  #将结果进行整合
  cbdRst = NULL
  # 1. 数据分布情况和中位生存时间
  survFitTb = summary(kmRst[['分组中位生存时间']])[["table"]]
  n = survFitTb[, 'records'] 
  disRatio = sapply(n*100 / sum(n), pointsZeroFill, pointsWidth = 1)
  nCmb = sprintf('%s (%s)', n, disRatio)
  #中位生存时间
  medSurvTime = ifelse(is.na(survFitTb[, 'median']), 'Not reached', sprintf('%s (%s, %s)', sapply(survFitTb[, 'median'], pointsZeroFill), sapply(survFitTb[, '0.95LCL'], pointsZeroFill), sapply(survFitTb[, '0.95UCL'], pointsZeroFill)))
  # 2. 单因素cox得到的HR和p值
  sgRst = kmRst[['单因素cox得到的HR等指标']]
  Hrs = sprintf('%s (%s, %s)', sapply(exp(sgRst[, 'coef']), pointsZeroFill), sapply(exp(sgRst[, 'coef'] - 1.96*sgRst[, 'se(coef)']), pointsZeroFill), sapply(exp(sgRst[, 'coef'] + 1.96*sgRst[, 'se(coef)']), pointsZeroFill))
  p = sapply(sgRst[, 'Pr(>|z|)'], singleP4Show)
  # 3. 指定时间点的生存概率
  tmSProb = kmRst[['指定时间点的生存概率']]
  nm = sprintf('%s#%s', tmSProb[, '指定时间点'], tmSProb[, '分层标识'])
  pbs = sprintf('%s (%s, %s)', sapply(as.numeric(tmSProb[, '生存概率'])*100, pointsZeroFill), sapply(as.numeric(tmSProb[, '生存概率_95ciLow'])*100, pointsZeroFill), sapply(as.numeric(tmSProb[, '生存概率_95ciUpper'])*100, pointsZeroFill))
  
  cbdRst = rbind(cbdRst
                  ,cbind(c('数据分布'), rownames(survFitTb), nCmb)
                  ,cbind(c('中位生存时间'), rownames(survFitTb), as.vector(medSurvTime))
                  ,cbind(c('HR'), rownames(sgRst), Hrs)
                  ,cbind(c('P'), rownames(sgRst), as.vector(p))
                  ,cbind(c('2年生存率'), nm, pbs)
  )
  return(cbdRst)
}
