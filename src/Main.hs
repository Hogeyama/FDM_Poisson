
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Main where

import           PoissonFDM                             hiding(Analysis(..))
import qualified PoissonFDM                             as P (Analysis(..))
import           Data.List                              (unzip4,intersperse)
import           Data.Colour                            (opaque)
import           Data.Colour.Names                      (blue,green,red)
import           Data.Default                           (def)
import           Control.Monad                          (void)
import           Control.Lens.Operators                 ((.~),(&))
import           System.IO                              (IOMode(..), hPutStrLn, withFile)
import           Graphics.Rendering.Chart
import           Graphics.Rendering.Chart.Backend.Cairo
import           GHC.TypeLits                           (KnownNat, type (^))
import           Control.Monad (forM_)

-------------------------------------------------------------------------------
-- example
-------------------------------------------------------------------------------

actual_u, laplacian :: Double -> Double -> Double
actual_u x y  = sin(3*x*y)
laplacian x y = -9*(x*x+y*y)*sin(3*x*y)

example :: forall n. KnownNat n => PoissonDef n
example = buildPoissonDef actual_u laplacian

-------------------------------------------------------------------------------
-- Data
-------------------------------------------------------------------------------

result :: [(Double, Double, Double, Double)]
result = map f [
    g (example :: PoissonDef 10)
  , g (example :: PoissonDef 20)
  , g (example :: PoissonDef 30)
  , g (example :: PoissonDef 40)
  , g (example :: PoissonDef 50)
  , g (example :: PoissonDef 60)
  , g (example :: PoissonDef 70)
  ]
  where
    f :: P.Analysis -> (Double,Double,Double,Double)
    f (P.Analysis n x y r) = (fromIntegral n, x, y, fromIntegral r)
    g :: forall n. (KnownNat n, KnownNat (n^2)) => PoissonDef n -> P.Analysis
    g = flip analyze actual_u

ns, consistency, convergence, repetition, hs:: [Double]
(ns, consistency, convergence, repetition) = unzip4 result
hs = map (\n -> 1.0 / (n+1)) ns

-------------------------------------------------------------------------------
-- Chart
-------------------------------------------------------------------------------

doubleLogLayout :: Layout Double Double
doubleLogLayout = def
  & layout_x_axis.laxis_generate .~ autoScaledLogAxis logAxisParams
  & layout_y_axis.laxis_generate .~ autoScaledLogAxis logAxisParams
  where logAxisParams = LogAxisParams $ map show

consistencyChart :: Renderable ()
consistencyChart = toRenderable layout
  where
    layout = doubleLogLayout
      & layout_title .~ "Consistensy"
      & layout_plots .~ [consistP,consistL,linear,quad]
      & layout_x_axis.laxis_title .~ "h"
      & layout_y_axis.laxis_title .~ ""
    consistencyData = zip hs consistency
    consistP = toPlot $ def
      & plot_points_values          .~ consistencyData
      & plot_points_style           .~ filledPolygon 5.0 4 False (opaque green)
      & plot_points_title           .~ "|b-Au*|_max"
    consistL = toPlot $ def
      & plot_lines_values           .~ [ consistencyData ]
      & plot_lines_style.line_color .~ opaque green
      & plot_lines_style.line_width .~ 3
      & plot_lines_title            .~ "|b-Au*|_max"
    linear = toPlot $ def
      & plot_lines_values           .~ [[ (x, 2*x) | x <- [0.01,0.10] ]]
      & plot_lines_style.line_color .~ opaque blue
      & plot_lines_style.line_width .~ 3
      & plot_lines_title            .~ "kh"
    quad = toPlot $ def
      & plot_lines_values           .~ [[ (x, 2*x*x) | x <- [0.01,0.10] ]]
      & plot_lines_style.line_color .~ opaque red
      & plot_lines_style.line_width .~ 3
      & plot_lines_title            .~ "kh^2"

convergenceChart :: Renderable ()
convergenceChart = toRenderable layout
  where
    layout = doubleLogLayout
      & layout_title .~ "Convergence"
      & layout_plots .~ [converP,converL,linear,quad]
      & layout_x_axis.laxis_title .~ "n"
      & layout_y_axis.laxis_title .~ ""
    convergenceData = zip hs convergence
    converP = toPlot $ def
      & plot_points_values          .~ convergenceData
      & plot_points_style           .~ filledPolygon 5.0 4 False (opaque green)
      & plot_points_title           .~ "|u-u*|_max"
    converL = toPlot $ def
      & plot_lines_values           .~ [ convergenceData ]
      & plot_lines_style.line_color .~ opaque green
      & plot_lines_style.line_width .~ 3
      & plot_lines_title            .~ "|u-u*|_max"
    linear = toPlot $ def
      & plot_lines_values           .~ [[ (x, 0.1*x) | x <- [0.01,0.02..0.10] ]]
      & plot_lines_style.line_color .~ opaque blue
      & plot_lines_style.line_width .~ 3
      & plot_lines_title            .~ "kh"
    quad = toPlot $ def
      & plot_lines_values           .~ [[ (x, 0.1*x*x) | x <- [0.01,0.10] ]]
      & plot_lines_style.line_color .~ opaque red
      & plot_lines_style.line_width .~ 3
      & plot_lines_title            .~ "kh^2"

repetitionChart :: Renderable ()
repetitionChart = toRenderable layout
  where
    layout = def
      & layout_title .~ "Repetition"
      & layout_plots .~ [repeatP,repeatL]
      & layout_x_axis.laxis_title .~ "n"
      & layout_y_axis.laxis_title .~ ""
    repetitionData = zip ns repetition
    repeatP = toPlot $ def
      & plot_points_values          .~ repetitionData
      & plot_points_style           .~ filledPolygon 5.0 4 False (opaque green)
      & plot_points_title           .~ "repetition"
    repeatL = toPlot $ def
      & plot_lines_values           .~ [ repetitionData ]
      & plot_lines_style.line_color .~ opaque green
      & plot_lines_style.line_width .~ 3
      & plot_lines_title            .~ "repetition"

-------------------------------------------------------------------------------
-- 
-------------------------------------------------------------------------------

main :: IO ()
main = do
  -- CSV for "plot3d.py"
  withFile "plot/plot.csv" WriteMode $ \h ->
    forM_ (answer $ solvePoisson (example :: PoissonDef 50)) $ \(x,y,z) ->
      hPutStrLn h $ unwords $ intersperse "," $ map show [x,y,z,actual_u x y]

  -- PNGs
  void $ renderableToFile def "./plot/consistency.png" consistencyChart
  void $ renderableToFile def "./plot/convergence.png" convergenceChart
  void $ renderableToFile def "./plot/repetition.png"  repetitionChart

