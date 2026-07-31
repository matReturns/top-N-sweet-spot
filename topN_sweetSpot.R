



# ============================================================
# Multi-Universe Top-N Concentration Study
# Core 5 vs Broader AAA 10 Universe
# ============================================================

library(quantmod)
library(PerformanceAnalytics)
library(xts)
library(zoo)

stratStats <- function(rets, digits = 4) {
  
  rets <- na.omit(rets)
  annual <- apply.yearly(rets, Return.cumulative)
  
  stats <- rbind(
    Return.annualized(rets),
    StdDev.annualized(rets),
    SharpeRatio.annualized(rets, Rf = 0),
    maxDrawdown(rets),
    CalmarRatio(rets),
    apply(annual, 2, min, na.rm = TRUE)
  )
  
  rownames(stats) <- c(
    "Annualized Return",
    "Annualized Volatility",
    "Annualized Sharpe",
    "Worst Drawdown",
    "Calmar Ratio",
    "Worst Calendar Year"
  )
  
  return(round(stats, digits))
}




run_multi_universe_concentration_test <- function(
    dataStartDate = "2006-01-01",
    analysisStartDate = "2007-07-01",
    endDate = Sys.Date(),
    universes = list(
      Core5 = c("SPY", "EFA", "EEM", "VNQ", "DBC"),
      AAA10 = c("SPY", "VGK", "EWJ", "EEM", "VNQ", "RWX", "IEF", "TLT", "DBC", "GLD")
    ),
    concentrationFractions = c(0.20, 0.40, 0.60, 0.80, 1.00),
    assetMomWindows = c(63, 126),
    canaryConfigs = list(
      EEM_BND_W21_63 = list(
        canaryAssets = c("EEM", "BND"),
        canaryMomWindows = c(21, 63)
      ),
      VWO_BND_W20_60 = list(
        canaryAssets = c("VWO", "BND"),
        canaryMomWindows = c(20, 60)
      )
    ),
    riskOffMode = c("conditional_defensive", "defensive", "cash"),
    defensiveAsset = "IEF",
    defensiveMomWindows = c(20, 60),
    cashAsset = NULL,
    cashReturn = 0,
    rebalanceOn = "months",
    reportYear = "2022",
    verbose = TRUE) {
  
  riskOffMode <- match.arg(riskOffMode)
  
  # -----------------------------
  # Helper functions
  # -----------------------------
  cumulative_return <- function(x) {
    prod(1 + as.numeric(x), na.rm = FALSE) - 1
  }
  
  momentum_score <- function(R, symbol, endRow, windows) {
    
    scores <- sapply(windows, function(w) {
      
      startRow <- endRow - w + 1
      
      if (startRow < 1) {
        return(NA_real_)
      }
      
      cumulative_return(R[startRow:endRow, symbol])
    })
    
    sum(scores, na.rm = FALSE)
  }
  
  get_year_return <- function(rets, year = "2022") {
    
    yearRange <- paste0(year, "-01-01/", year, "-12-31")
    yr <- rets[yearRange]
    
    if (NROW(yr) == 0) {
      return(NA_real_)
    }
    
    as.numeric(Return.cumulative(yr))
  }
  
  round_numeric_df <- function(df, digits = 4) {
    
    numericCols <- sapply(df, is.numeric)
    df[numericCols] <- lapply(df[numericCols], round, digits)
    
    return(df)
  }
  
  make_selection_label <- function(frac) {
    paste0("Top", round(frac * 100), "Pct")
  }
  
  # -----------------------------
  # Symbols needed
  # -----------------------------
  allUniverseSymbols <- unique(unlist(universes))
  canarySymbols <- unique(unlist(lapply(canaryConfigs, function(x) x$canaryAssets)))
  
  downloadSymbols <- unique(c(
    allUniverseSymbols,
    canarySymbols,
    defensiveAsset,
    cashAsset
  ))
  
  downloadSymbols <- downloadSymbols[!is.na(downloadSymbols)]
  
  if (is.null(cashAsset)) {
    cashName <- "CASH"
  } else {
    cashName <- cashAsset
  }
  
  # -----------------------------
  # Download adjusted prices
  # -----------------------------
  priceList <- list()
  
  for (sym in downloadSymbols) {
    
    if (verbose) {
      message("Downloading: ", sym)
    }
    
    tmp <- getSymbols(
      sym,
      from = dataStartDate,
      to = endDate,
      auto.assign = FALSE,
      warnings = FALSE
    )
    
    px <- Ad(tmp)
    colnames(px) <- sym
    priceList[[sym]] <- px
  }
  
  prices <- do.call(merge, priceList)
  prices <- na.omit(prices)
  
  assetReturns <- Return.calculate(prices)
  assetReturns <- na.omit(assetReturns)
  
  if (is.null(cashAsset)) {
    
    cashReturns <- xts(
      rep(cashReturn / 252, NROW(assetReturns)),
      order.by = index(assetReturns)
    )
    
    colnames(cashReturns) <- "CASH"
    
    allReturns <- merge(assetReturns, cashReturns)
    
  } else {
    
    allReturns <- assetReturns
  }
  
  allReturns <- na.omit(allReturns)
  
  if (verbose) {
    message("Return data begins: ", as.character(first(index(allReturns))))
    message("Return data ends:   ", as.character(last(index(allReturns))))
  }
  
  # -----------------------------
  # Rebalance dates
  # -----------------------------
  allCanaryWindows <- unique(unlist(lapply(canaryConfigs, function(x) x$canaryMomWindows)))
  
  maxLookback <- max(
    assetMomWindows,
    allCanaryWindows,
    defensiveMomWindows
  )
  
  rebalanceRows <- endpoints(allReturns, on = rebalanceOn)
  rebalanceRows <- rebalanceRows[rebalanceRows > 0]
  rebalanceRows <- rebalanceRows[rebalanceRows <= NROW(allReturns)]
  signalRows <- rebalanceRows[rebalanceRows > maxLookback]
  
  # ============================================================
  # Equal Weight Benchmark
  # ============================================================
  run_equal_weight <- function(universeName, investableAssets) {
    
    rets <- Return.portfolio(
      R = allReturns[, investableAssets],
      weights = rep(1 / length(investableAssets), length(investableAssets)),
      rebalance_on = rebalanceOn
    )
    
    colnames(rets) <- paste0(universeName, "_EqualWeight")
    
    rets <- rets[paste0(analysisStartDate, "/")]
    rets <- na.omit(rets)
    
    return(rets)
  }
  
  # ============================================================
  # Unfiltered Momentum
  # ============================================================
  run_unfiltered_momentum <- function(universeName,
                                      investableAssets,
                                      nAssets,
                                      selectionLabel,
                                      targetFraction) {
    
    strategyName <- paste(
      universeName,
      selectionLabel,
      "Momentum",
      sep = "_"
    )
    
    weightsMat <- matrix(
      NA_real_,
      nrow = NROW(allReturns),
      ncol = length(investableAssets)
    )
    
    colnames(weightsMat) <- investableAssets
    
    for (signalRow in signalRows) {
      
      assetScores <- sapply(investableAssets, function(sym) {
        momentum_score(allReturns, sym, signalRow, assetMomWindows)
      })
      
      rankedAssets <- names(sort(assetScores, decreasing = TRUE))
      selectedAssets <- rankedAssets[1:nAssets]
      
      w <- rep(0, length(investableAssets))
      names(w) <- investableAssets
      w[selectedAssets] <- 1 / nAssets
      
      weightsMat[signalRow, ] <- w[investableAssets]
    }
    
    weights <- xts(weightsMat, order.by = index(allReturns))
    weights <- na.locf(weights, na.rm = FALSE)
    weightsLag <- lag(weights, k = 1)
    
    rets <- xts(
      rowSums(weightsLag * allReturns[, investableAssets], na.rm = FALSE),
      order.by = index(allReturns)
    )
    
    colnames(rets) <- strategyName
    
    retsLive <- rets[paste0(analysisStartDate, "/")]
    retsLive <- na.omit(retsLive)
    
    liveWeights <- weightsLag[index(retsLive), investableAssets]
    
    weightChanges <- abs(liveWeights - lag(liveWeights, k = 1))
    
    dailyTurnover <- xts(
      rowSums(weightChanges, na.rm = TRUE),
      order.by = index(liveWeights)
    )
    
    annualTurnover <- apply.yearly(dailyTurnover, sum, na.rm = TRUE)
    
    turnoverSummary <- data.frame(
      Strategy = strategyName,
      Universe = universeName,
      StrategyType = "Unfiltered Momentum",
      SelectionLabel = selectionLabel,
      TargetFraction = targetFraction,
      nAssets = nAssets,
      UniverseSize = length(investableAssets),
      AverageAnnualTurnover = mean(as.numeric(annualTurnover), na.rm = TRUE),
      MaximumAnnualTurnover = max(as.numeric(annualTurnover), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    return(list(
      returns = retsLive,
      weights = weights,
      weightsLag = weightsLag,
      turnover = turnoverSummary
    ))
  }
  
  # ============================================================
  # Canary Momentum
  # ============================================================
  run_canary_momentum <- function(universeName,
                                  investableAssets,
                                  nAssets,
                                  selectionLabel,
                                  targetFraction,
                                  configName,
                                  canaryAssets,
                                  canaryMomWindows) {
    
    strategyName <- paste(
      universeName,
      selectionLabel,
      configName,
      sep = "_"
    )
    
    tradeAssets <- unique(c(
      investableAssets,
      defensiveAsset,
      cashName
    ))
    
    weightsMat <- matrix(
      NA_real_,
      nrow = NROW(allReturns),
      ncol = length(tradeAssets)
    )
    
    colnames(weightsMat) <- tradeAssets
    
    riskBudgetVec <- rep(NA_real_, NROW(allReturns))
    
    for (signalRow in signalRows) {
      
      # -----------------------------
      # Rank investable assets
      # -----------------------------
      assetScores <- sapply(investableAssets, function(sym) {
        momentum_score(allReturns, sym, signalRow, assetMomWindows)
      })
      
      rankedAssets <- names(sort(assetScores, decreasing = TRUE))
      selectedAssets <- rankedAssets[1:nAssets]
      
      # -----------------------------
      # Canary risk budget
      # -----------------------------
      canaryScores <- sapply(canaryAssets, function(sym) {
        momentum_score(allReturns, sym, signalRow, canaryMomWindows)
      })
      
      riskBudget <- mean(canaryScores > 0)
      riskOffWeight <- 1 - riskBudget
      
      # -----------------------------
      # Start weights
      # -----------------------------
      w <- rep(0, length(tradeAssets))
      names(w) <- tradeAssets
      
      # Risk-on allocation
      w[selectedAssets] <- riskBudget / nAssets
      
      # Risk-off allocation
      if (riskOffWeight > 0) {
        
        if (riskOffMode == "cash") {
          
          w[cashName] <- w[cashName] + riskOffWeight
          
        } else if (riskOffMode == "defensive") {
          
          w[defensiveAsset] <- w[defensiveAsset] + riskOffWeight
          
        } else if (riskOffMode == "conditional_defensive") {
          
          defensiveScore <- momentum_score(
            allReturns,
            defensiveAsset,
            signalRow,
            defensiveMomWindows
          )
          
          if (!is.na(defensiveScore) && defensiveScore > 0) {
            w[defensiveAsset] <- w[defensiveAsset] + riskOffWeight
          } else {
            w[cashName] <- w[cashName] + riskOffWeight
          }
        }
      }
      
      weightsMat[signalRow, ] <- w[tradeAssets]
      riskBudgetVec[signalRow] <- riskBudget
    }
    
    weights <- xts(weightsMat, order.by = index(allReturns))
    weights <- na.locf(weights, na.rm = FALSE)
    weightsLag <- lag(weights, k = 1)
    
    rets <- xts(
      rowSums(weightsLag * allReturns[, tradeAssets], na.rm = FALSE),
      order.by = index(allReturns)
    )
    
    colnames(rets) <- strategyName
    
    retsLive <- rets[paste0(analysisStartDate, "/")]
    retsLive <- na.omit(retsLive)
    
    liveWeights <- weightsLag[index(retsLive), tradeAssets]
    
    riskBudget <- xts(riskBudgetVec, order.by = index(allReturns))
    riskBudget <- na.locf(riskBudget, na.rm = FALSE)
    riskBudgetLag <- lag(riskBudget, k = 1)
    liveRiskBudget <- riskBudgetLag[index(retsLive)]
    
    # -----------------------------
    # Exposure
    # -----------------------------
    exposure <- data.frame(
      Strategy = strategyName,
      Universe = universeName,
      SelectionLabel = selectionLabel,
      TargetFraction = targetFraction,
      nAssets = nAssets,
      UniverseSize = length(investableAssets),
      CanaryConfig = configName,
      avgRiskBudget = mean(as.numeric(liveRiskBudget), na.rm = TRUE),
      pctFullRiskOn = mean(as.numeric(liveRiskBudget) == 1, na.rm = TRUE),
      pctHalfRiskOn = mean(as.numeric(liveRiskBudget) == 0.5, na.rm = TRUE),
      pctRiskOff = mean(as.numeric(liveRiskBudget) == 0, na.rm = TRUE),
      avgRiskAssetWeight = mean(rowSums(liveWeights[, investableAssets], na.rm = TRUE), na.rm = TRUE),
      avgDefensiveWeight = if (defensiveAsset %in% colnames(liveWeights)) {
        mean(as.numeric(liveWeights[, defensiveAsset]), na.rm = TRUE)
      } else {
        NA_real_
      },
      avgCashWeight = mean(as.numeric(liveWeights[, cashName]), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    # -----------------------------
    # Turnover
    # -----------------------------
    weightChanges <- abs(liveWeights - lag(liveWeights, k = 1))
    
    dailyTurnover <- xts(
      rowSums(weightChanges, na.rm = TRUE),
      order.by = index(liveWeights)
    )
    
    annualTurnover <- apply.yearly(dailyTurnover, sum, na.rm = TRUE)
    
    turnoverSummary <- data.frame(
      Strategy = strategyName,
      Universe = universeName,
      StrategyType = "Canary Momentum",
      SelectionLabel = selectionLabel,
      TargetFraction = targetFraction,
      nAssets = nAssets,
      UniverseSize = length(investableAssets),
      CanaryConfig = configName,
      AverageAnnualTurnover = mean(as.numeric(annualTurnover), na.rm = TRUE),
      MaximumAnnualTurnover = max(as.numeric(annualTurnover), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    return(list(
      returns = retsLive,
      weights = weights,
      weightsLag = weightsLag,
      riskBudget = riskBudget,
      exposure = exposure,
      turnover = turnoverSummary
    ))
  }
  
  # ============================================================
  # Run all universes and selection fractions
  # ============================================================
  equalWeightList <- list()
  unfilteredResults <- list()
  unfilteredReturnsList <- list()
  unfilteredTurnoverList <- list()
  
  canaryResults <- list()
  canaryReturnsList <- list()
  exposureList <- list()
  canaryTurnoverList <- list()
  
  for (universeName in names(universes)) {
    
    investableAssets <- universes[[universeName]]
    universeSize <- length(investableAssets)
    
    if (verbose) {
      message("Running universe: ", universeName)
    }
    
    # Equal weight benchmark
    ew <- run_equal_weight(universeName, investableAssets)
    equalWeightList[[colnames(ew)]] <- ew
    
    for (frac in concentrationFractions) {
      
      selectionLabel <- make_selection_label(frac)
      nAssets <- round(frac * universeSize)
      nAssets <- max(1, min(universeSize, nAssets))
      
      if (verbose) {
        message("  Selection: ", selectionLabel, " | Top ", nAssets, " of ", universeSize)
      }
      
      # Unfiltered momentum
      unf <- run_unfiltered_momentum(
        universeName = universeName,
        investableAssets = investableAssets,
        nAssets = nAssets,
        selectionLabel = selectionLabel,
        targetFraction = frac
      )
      
      unfilteredResults[[colnames(unf$returns)]] <- unf
      unfilteredReturnsList[[colnames(unf$returns)]] <- unf$returns
      unfilteredTurnoverList[[colnames(unf$returns)]] <- unf$turnover
      
      # Canary versions
      for (configName in names(canaryConfigs)) {
        
        res <- run_canary_momentum(
          universeName = universeName,
          investableAssets = investableAssets,
          nAssets = nAssets,
          selectionLabel = selectionLabel,
          targetFraction = frac,
          configName = configName,
          canaryAssets = canaryConfigs[[configName]]$canaryAssets,
          canaryMomWindows = canaryConfigs[[configName]]$canaryMomWindows
        )
        
        canaryResults[[colnames(res$returns)]] <- res
        canaryReturnsList[[colnames(res$returns)]] <- res$returns
        exposureList[[colnames(res$returns)]] <- res$exposure
        canaryTurnoverList[[colnames(res$returns)]] <- res$turnover
      }
    }
  }
  
  equalWeightReturns <- do.call(merge, equalWeightList)
  unfilteredReturns <- do.call(merge, unfilteredReturnsList)
  canaryReturns <- do.call(merge, canaryReturnsList)
  
  exposureStats <- do.call(rbind, exposureList)
  unfilteredTurnover <- do.call(rbind, unfilteredTurnoverList)
  canaryTurnover <- do.call(rbind, canaryTurnoverList)
  
  rownames(exposureStats) <- NULL
  rownames(unfilteredTurnover) <- NULL
  rownames(canaryTurnover) <- NULL
  
  # -----------------------------
  # Combine returns
  # -----------------------------
  mainReturns <- merge(
    equalWeightReturns,
    unfilteredReturns,
    canaryReturns
  )
  
  mainReturns <- mainReturns[paste0(analysisStartDate, "/")]
  mainReturns <- na.omit(mainReturns)
  
  canaryReturns <- canaryReturns[index(mainReturns)]
  canaryReturns <- na.omit(canaryReturns)
  
  unfilteredReturns <- unfilteredReturns[index(mainReturns)]
  unfilteredReturns <- na.omit(unfilteredReturns)
  
  # -----------------------------
  # Summary
  # -----------------------------
  summary <- stratStats(mainReturns)
  
  summaryTable <- data.frame(
    Strategy = colnames(mainReturns),
    AnnualizedReturn = as.numeric(summary["Annualized Return", ]),
    AnnualizedVolatility = as.numeric(summary["Annualized Volatility", ]),
    Sharpe = as.numeric(summary["Annualized Sharpe", ]),
    WorstDrawdown = as.numeric(summary["Worst Drawdown", ]),
    Calmar = as.numeric(summary["Calmar Ratio", ]),
    WorstCalendarYear = as.numeric(summary["Worst Calendar Year", ]),
    stringsAsFactors = FALSE
  )
  
  yearColName <- paste0("Return", reportYear)
  
  summaryTable[[yearColName]] <- sapply(
    summaryTable$Strategy,
    function(s) get_year_return(mainReturns[, s], reportYear)
  )
  
  # -----------------------------
  # Add metadata
  # -----------------------------
  summaryTable$Universe <- NA_character_
  summaryTable$StrategyType <- NA_character_
  summaryTable$SelectionLabel <- NA_character_
  summaryTable$TargetFraction <- NA_real_
  summaryTable$nAssets <- NA_integer_
  summaryTable$UniverseSize <- NA_integer_
  summaryTable$CanaryConfig <- NA_character_
  
  # Equal weight metadata
  for (universeName in names(universes)) {
    
    ewName <- paste0(universeName, "_EqualWeight")
    
    idx <- summaryTable$Strategy == ewName
    
    summaryTable$Universe[idx] <- universeName
    summaryTable$StrategyType[idx] <- "Equal Weight"
    summaryTable$SelectionLabel[idx] <- "All"
    summaryTable$TargetFraction[idx] <- 1
    summaryTable$nAssets[idx] <- length(universes[[universeName]])
    summaryTable$UniverseSize[idx] <- length(universes[[universeName]])
  }
  
  # Unfiltered metadata
  for (nm in names(unfilteredResults)) {
    
    idx <- summaryTable$Strategy == nm
    meta <- unfilteredResults[[nm]]$turnover
    
    summaryTable$Universe[idx] <- meta$Universe
    summaryTable$StrategyType[idx] <- "Unfiltered Momentum"
    summaryTable$SelectionLabel[idx] <- meta$SelectionLabel
    summaryTable$TargetFraction[idx] <- meta$TargetFraction
    summaryTable$nAssets[idx] <- meta$nAssets
    summaryTable$UniverseSize[idx] <- meta$UniverseSize
  }
  
  # Canary metadata
  for (nm in names(canaryResults)) {
    
    idx <- summaryTable$Strategy == nm
    meta <- canaryResults[[nm]]$turnover
    
    summaryTable$Universe[idx] <- meta$Universe
    summaryTable$StrategyType[idx] <- "Canary Momentum"
    summaryTable$SelectionLabel[idx] <- meta$SelectionLabel
    summaryTable$TargetFraction[idx] <- meta$TargetFraction
    summaryTable$nAssets[idx] <- meta$nAssets
    summaryTable$UniverseSize[idx] <- meta$UniverseSize
    summaryTable$CanaryConfig[idx] <- meta$CanaryConfig
  }
  
  summaryTable <- summaryTable[, c(
    "Strategy",
    "Universe",
    "StrategyType",
    "SelectionLabel",
    "TargetFraction",
    "nAssets",
    "UniverseSize",
    "CanaryConfig",
    "AnnualizedReturn",
    "AnnualizedVolatility",
    "Sharpe",
    "WorstDrawdown",
    "Calmar",
    "WorstCalendarYear",
    yearColName
  )]
  
  # -----------------------------
  # Rankings
  # -----------------------------
  rankedBySharpe <- summaryTable[order(-summaryTable$Sharpe), ]
  rankedByCalmar <- summaryTable[order(-summaryTable$Calmar), ]
  rankedByDrawdown <- summaryTable[order(summaryTable$WorstDrawdown), ]
  rankedByReturn <- summaryTable[order(-summaryTable$AnnualizedReturn), ]
  
  rownames(rankedBySharpe) <- NULL
  rownames(rankedByCalmar) <- NULL
  rownames(rankedByDrawdown) <- NULL
  rownames(rankedByReturn) <- NULL
  
  annualReturns <- apply.yearly(mainReturns, Return.cumulative)
  
  # -----------------------------
  # Aggregates
  # -----------------------------
  canaryOnly <- subset(summaryTable, StrategyType == "Canary Momentum")
  
  selectionSummary <- aggregate(
    cbind(
      AnnualizedReturn,
      AnnualizedVolatility,
      Sharpe,
      WorstDrawdown,
      Calmar,
      WorstCalendarYear
    ) ~ Universe + SelectionLabel + TargetFraction + nAssets + UniverseSize,
    data = canaryOnly,
    FUN = mean
  )
  
  selectionSummary <- selectionSummary[order(
    selectionSummary$Universe,
    selectionSummary$TargetFraction
  ), ]
  
  rownames(selectionSummary) <- NULL
  
  universeSummary <- aggregate(
    cbind(
      AnnualizedReturn,
      AnnualizedVolatility,
      Sharpe,
      WorstDrawdown,
      Calmar,
      WorstCalendarYear
    ) ~ Universe,
    data = canaryOnly,
    FUN = mean
  )
  
  universeSummary <- universeSummary[order(-universeSummary$Sharpe), ]
  rownames(universeSummary) <- NULL
  
  configSummary <- aggregate(
    cbind(
      AnnualizedReturn,
      AnnualizedVolatility,
      Sharpe,
      WorstDrawdown,
      Calmar,
      WorstCalendarYear
    ) ~ Universe + CanaryConfig,
    data = canaryOnly,
    FUN = mean
  )
  
  configSummary <- configSummary[order(configSummary$Universe, -configSummary$Sharpe), ]
  rownames(configSummary) <- NULL
  
  # Matrices based on averaged canary results
  sharpeMatrix <- xtabs(
    Sharpe ~ SelectionLabel + Universe,
    data = selectionSummary
  )
  
  calmarMatrix <- xtabs(
    Calmar ~ SelectionLabel + Universe,
    data = selectionSummary
  )
  
  drawdownMatrix <- xtabs(
    WorstDrawdown ~ SelectionLabel + Universe,
    data = selectionSummary
  )
  
  return(list(
    mainReturns = mainReturns,
    canaryReturns = canaryReturns,
    unfilteredReturns = unfilteredReturns,
    equalWeightReturns = equalWeightReturns,
    summary = summary,
    summaryTable = round_numeric_df(summaryTable, 4),
    rankedBySharpe = round_numeric_df(rankedBySharpe, 4),
    rankedByCalmar = round_numeric_df(rankedByCalmar, 4),
    rankedByDrawdown = round_numeric_df(rankedByDrawdown, 4),
    rankedByReturn = round_numeric_df(rankedByReturn, 4),
    selectionSummary = round_numeric_df(selectionSummary, 4),
    universeSummary = round_numeric_df(universeSummary, 4),
    configSummary = round_numeric_df(configSummary, 4),
    exposureStats = round_numeric_df(exposureStats, 4),
    unfilteredTurnover = round_numeric_df(unfilteredTurnover, 4),
    canaryTurnover = round_numeric_df(canaryTurnover, 4),
    annualReturns = round(annualReturns, 4),
    sharpeMatrix = round(sharpeMatrix, 4),
    calmarMatrix = round(calmarMatrix, 4),
    drawdownMatrix = round(drawdownMatrix, 4),
    canaryResults = canaryResults,
    unfilteredResults = unfilteredResults,
    prices = prices,
    allReturns = allReturns,
    settings = list(
      dataStartDate = dataStartDate,
      analysisStartDate = as.character(first(index(mainReturns))),
      endDate = as.character(last(index(mainReturns))),
      universes = universes,
      concentrationFractions = concentrationFractions,
      assetMomWindows = assetMomWindows,
      canaryConfigs = canaryConfigs,
      riskOffMode = riskOffMode,
      defensiveAsset = defensiveAsset,
      defensiveMomWindows = defensiveMomWindows,
      cashAsset = cashAsset,
      cashReturn = cashReturn,
      rebalanceOn = rebalanceOn,
      reportYear = reportYear
    )
  ))
}




multiUniverseTest <- run_multi_universe_concentration_test(
  dataStartDate = "2006-01-01",
  analysisStartDate = "2007-07-01",
  endDate = "2026-07-09",
  universes = list(
    Core5 = c("SPY", "EFA", "EEM", "VNQ", "DBC"),
    AAA10 = c("SPY", "VGK", "EWJ", "EEM", "VNQ", "RWX", "IEF", "TLT", "DBC", "GLD")
  ),
  concentrationFractions = c(0.20, 0.40, 0.60, 0.80, 1.00),
  assetMomWindows = c(63, 126),
  canaryConfigs = list(
    EEM_BND_W21_63 = list(
      canaryAssets = c("EEM", "BND"),
      canaryMomWindows = c(21, 63)
    ),
    VWO_BND_W20_60 = list(
      canaryAssets = c("VWO", "BND"),
      canaryMomWindows = c(20, 60)
    )
  ),
  riskOffMode = "conditional_defensive",
  defensiveAsset = "IEF",
  defensiveMomWindows = c(20, 60),
  cashAsset = NULL,
  cashReturn = 0,
  rebalanceOn = "months",
  reportYear = "2022",
  verbose = TRUE
)




multiUniverseTest$summaryTable
multiUniverseTest$rankedBySharpe
multiUniverseTest$rankedByCalmar
multiUniverseTest$rankedByDrawdown
multiUniverseTest$rankedByReturn

multiUniverseTest$selectionSummary
multiUniverseTest$universeSummary
multiUniverseTest$configSummary

multiUniverseTest$sharpeMatrix
multiUniverseTest$calmarMatrix
multiUniverseTest$drawdownMatrix

multiUniverseTest$exposureStats
multiUniverseTest$canaryTurnover
multiUniverseTest$unfilteredTurnover
multiUniverseTest$annualReturns
multiUniverseTest$settings




top6SharpeNames <- multiUniverseTest$rankedBySharpe$Strategy[
  multiUniverseTest$rankedBySharpe$StrategyType == "Canary Momentum"
][1:6]

charts.PerformanceSummary(
  multiUniverseTest$mainReturns[, top6SharpeNames],
  main = "Core 5 vs AAA 10: Best Results by Sharpe",
  wealth.index = T
)




top6CalmarNames <- multiUniverseTest$rankedByCalmar$Strategy[
  multiUniverseTest$rankedByCalmar$StrategyType == "Canary Momentum"
][1:6]

charts.PerformanceSummary(
  multiUniverseTest$mainReturns[, top6CalmarNames],
  main = "Core 5 vs AAA 10: Best Results by Calmar"
)




barplot(
  t(multiUniverseTest$sharpeMatrix),
  beside = TRUE,
  legend.text = TRUE,
  args.legend = list(x = "topright"),
  main = "Average Sharpe by Concentration Level",
  ylab = "Average Sharpe Across Canary Configurations",
  xlab = "Concentration Level"
)




barplot(
  t(multiUniverseTest$calmarMatrix),
  beside = TRUE,
  legend.text = TRUE,
  args.legend = list(x = "topright"),
  main = "Average Calmar by Concentration Level",
  ylab = "Average Calmar Across Canary Configurations",
  xlab = "Concentration Level"
)




barplot(
  t(multiUniverseTest$drawdownMatrix),
  beside = TRUE,
  legend.text = TRUE,
  args.legend = list(x = "topright"),
  main = "Average Worst Drawdown by Concentration Level",
  ylab = "Average Worst Drawdown",
  xlab = "Concentration Level"
)




multiUniverseTest$selectionSummary
multiUniverseTest$universeSummary
multiUniverseTest$rankedBySharpe
multiUniverseTest$rankedByCalmar
multiUniverseTest$sharpeMatrix
multiUniverseTest$calmarMatrix
multiUniverseTest$canaryTurnover







































