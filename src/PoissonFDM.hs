{-# LANGUAGE GADTs #-}
{-# LANGUAGE Strict #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE StrictData #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -Wall #-}

module PoissonFDM (
    PoissonDef(..)
  , PoissonResult(..)
  , Analysis(..)
  , SMat
  , SVec
  , printSMat
  , analyze
  , solvePoisson
  , cgMethod
  , buildPoissonDef
  , fromListVec
  , fromList
  , fromNat'
  ) where

import Data.Sized.Fin           (Fin, Sing, sing, universe, fromNat)
import Data.Sized.Matrix        (show2D)
import Data.Sized.Sparse.Matrix as SP
import GHC.TypeLits             (Nat, KnownNat, type (^))
import Data.List                ( unfoldr)

-------------------------------------------------------------------------------
-- Support function
-------------------------------------------------------------------------------

-- | 型レベル自然数をNumに変換
fromNat' :: forall f n a. (KnownNat n, Num a) => f n -> a
fromNat' _ = fromIntegral $ fromNat (sing :: Sing n)

-------------------------------------------------------------------------------
-- Sparse Matrix & Utility
-------------------------------------------------------------------------------

-- | n by m sparse matrix
--   Fin n は0以上n未満の自然数を表す
type SMat (n::Nat) (m::Nat) = SpMatrix (Fin n, Fin m) Double

-- | Sparse Vector
--   ベクトルは疎でないが，普通の行列と疎行列の積はライブラリにないのでともに疎なものを扱う
type SVec (n::Nat) = SMat n 1

--------------
-- 行列演算 --
--------------

-- | 行列積 A B
(<>) :: (KnownNat n, KnownNat k, KnownNat k', KnownNat m, k~k')
     => SMat n k -> SMat k' m -> SMat n m
(<>) = SP.mm
infixr 6 <>

-- | 行列積 A^T B
--(^<>) :: (KnownNat n, KnownNat k, KnownNat k', KnownNat m, k~k')
--      => SMat k n -> SMat k' m -> SMat n m
--a ^<> b = SP.mm (SP.transpose a) b
--infixr 6 ^<>

-- | 行列積 A B^T
--(<>^) :: (KnownNat n, KnownNat k, KnownNat k', KnownNat m, k~k')
--      => SMat n k -> SMat m k' -> SMat n m
--a <>^ b = SP.mm a (SP.transpose b)
--infixr 6 <>^

-- | 行列和
(.+.) :: (KnownNat n, KnownNat m) => SMat n m -> SMat n m -> SMat n m
(.+.) = SP.zipWith (+)
infixr 5 .+.

-- | 行列差
(.-.) :: (KnownNat n, KnownNat m) => SMat n m -> SMat n m -> SMat n m
(.-.) = SP.zipWith (-)
infixr 5 .-.

-- | 内積
(<.>) :: (KnownNat n) => SVec n -> SVec n -> Double
u <.> v = (transpose u <> v) ! (0 :: Int)

-- | u*Av
--   (u<#A#>v) で u*Avになる
(<#) :: KnownNat n => SVec n -> SVec n -> Double
(<#) = (<.>)
(#>) :: KnownNat n => SMat n n -> SVec n -> SVec n
(#>) = (<>)
infixr 4 <#, #>

-- | Frobenius norm (or Euclidean norm for Vector)
frobNorm :: (KnownNat n, KnownNat m) => SMat n m -> Double
frobNorm = sum . map (sq.snd) . snd . toAssocList
  where sq x = x*x

-- | max norm
maxNorm :: (KnownNat n, KnownNat m) => SMat n m -> Double
maxNorm = maximum . map (abs.snd) . snd . toAssocList

-- | zero matrix
zero :: (KnownNat n, KnownNat m) => SMat n m
zero = fromList []

---------------
-- from list --
---------------

fromList :: (KnownNat n, KnownNat m) => [((Fin n, Fin m), Double)] -> SMat n m
fromList = fromAssocList 0.0

fromListVec :: (KnownNat n) => [(Fin n, Double)] -> SVec n
fromListVec xs = fromList [ ((n,0),x) | (n,x) <- xs ]

-- | pretty print
--   密行列に変換するので遅いはず. デバグ用
printSMat :: (KnownNat n, KnownNat m) => SMat n m -> IO ()
printSMat = putStr . show2D . fill

-- | get element of vector
(!) :: (KnownNat n, Integral a) => SVec n -> a -> Double
(!) v i = getElem v (fromIntegral i, 0)

-- | scalar multiplication
(*^) :: (KnownNat n, KnownNat m) => Double -> SMat n m -> SMat n m
s *^ m = fmap (s*) m
infix 7 *^

-------------------------------------------------------------------------------
-- CG method
-------------------------------------------------------------------------------

-- (x,r,p) = (近似解, 残差, 勾配ベクトル)
type Triple (n :: Nat) = (SVec n, SVec n, SVec n)

cgMethod :: forall n. KnownNat n => SMat n n -> SVec n -> [Triple n]
cgMethod a b = (x0,r0,p0) : unfoldr step (x0,r0,p0)
  where
    eps  = 1e-12 :: Double

    x0, r0, p0 :: SVec n
    x0 = zero
    r0 = b .-. a<>x0
    p0 = r0

    step :: Triple n -> Maybe (Triple n, Triple n)
    step (x,r,p)
      | zansa < eps * frobNorm b = Nothing -- 相対残差eps以下
      | otherwise                = Just (triple,triple)
      where
        zansa  = frobNorm r
        alpha  = r<.>r / (p<#a#>p)
        beta   = r'<.>r' / r<.>r
        x'     = x .+. alpha *^ p
        r'     = r .-. alpha *^ a <> p
        p'     = r'.+. beta  *^ p
        triple = (x',r',p')

-------------------------------------------------------------------------------
-- Poisson
-------------------------------------------------------------------------------

data PoissonDef (n::Nat) = KnownNat n => PoissonDef {
    rangeX   :: (Double, Double)           -- ^ X-range
  , rangeY   :: (Double, Double)           -- ^ Y-range
  , f        :: Double -> Double -> Double
  , marginX0 :: SVec n                     -- ^ u(x_i,  y_-1)
  , margin0Y :: SVec n                     -- ^ u(x_-1, y_j )
  , marginXR :: SVec n                     -- ^ u(x_i,  y_n )
  , marginRY :: SVec n                     -- ^ u(x_n,  y_j )
  }

-- | i番目のx座標 (i=0~n-1)
ith_x :: forall n. PoissonDef n -> Fin n -> Double
ith_x PoissonDef{..} i = x_min + (x_max-x_min) * (fromIntegral i + 1.0) / (n + 1.0)
  where n = fromNat' (sing :: Sing n)
        (x_min,x_max) = rangeX

-- | i番目のy座標 (i=0~n-1)
ith_y :: forall n. PoissonDef n -> Fin n -> Double
ith_y PoissonDef{..} i = y_min + (y_max-y_min) * (fromIntegral i + 1.0) / (n + 1.0)
  where n = fromNat' (sing :: Sing n)
        (y_min,y_max) = rangeY

-- | (n*n) by (n*n) matrix
poissonMat :: forall n. (KnownNat n, KnownNat (n^2)) => PoissonDef n -> SMat (n^2) (n^2)
poissonMat _ = fromList xs
  where
    sn2 = (fromNat' (sing :: Sing n) + 1) ^ (2::Int) -- (n+1)^2
    xs = [ ((pos (i,j), pos (p,q)), x)
         | (i,j) <- universe
         , (p,q) <- neighbor (i,j)
         , let x | i==p && j==q =  4.0*sn2
                 | otherwise    = -1.0*sn2
         ]
    neighbor (i,j) = (i,j) : [ (i,q) | q <- nb j ] ++ [ (p,j) | p <- nb i ]
      where nb k | k==minBound = [k+1]
                 | k==maxBound = [k-1]
                 | otherwise   = [k-1, k+1]


-- | (n*n) vector
poissonVec :: forall n. (KnownNat n, KnownNat (n^2)) => PoissonDef n -> SVec (n^2)
poissonVec pp@PoissonDef{..} = fromListVec xs
  where
    sn2 = (fromNat' (sing :: Sing n) + 1) ^ (2 :: Int) -- (n+1)^2
    xs = [ (pos (i,j), b)
         | (i,j) <- universe
         , let b = - f (pp `ith_x` i) (pp `ith_y` j)
                   + sn2 * (if | i == minBound -> margin0Y ! j
                               | i == maxBound -> marginRY ! j
                               | otherwise -> 0.0)
                   + sn2 * (if | j == minBound -> marginX0 ! i
                               | j == maxBound -> marginXR ! i
                               | otherwise -> 0.0)
         ]

-- | u_{i,j}を一列に並べる
pos :: forall n. (KnownNat n, KnownNat (n^2)) => (Fin n,Fin n) -> Fin (n^2)
pos (i,j) = fromIntegral $ n*j'+i'
  where n  = fromNat' (sing :: Sing n) :: Int
        i' = fromIntegral i
        j' = fromIntegral j

----------------
-- CG法で解く --
----------------

data PoissonResult n where
  PoissonResult :: (KnownNat n, KnownNat (n^2)) => {
      triples :: [Triple (n^2)]
    , answer  :: [(Double,Double,Double)] -- (x,y,u)
    } -> PoissonResult n


solvePoisson :: forall n. (KnownNat n, KnownNat (n^2)) => PoissonDef n -> PoissonResult n
solvePoisson pp = PoissonResult{..}
  where
    a       = poissonMat pp
    b       = poissonVec pp
    triples = cgMethod a b
    answer  = [ (pp `ith_x` i, pp `ith_y` j, u!ix)
              | (i,j) <- universe
              , let ix = pos (i,j)
                    (u,_,_) = last triples
              ]

-------------------------------------------------------------------------------
-- Analysis
-------------------------------------------------------------------------------

data Analysis = Analysis {
    n           :: Int
  , consistency :: Double
  , convergence :: Double
  , repetition  :: Int
  }

analyze :: forall n. (KnownNat n, KnownNat (n^2)) => PoissonDef n -> (Double -> Double -> Double) -> Analysis
analyze pp exact_u = Analysis n x y (length triples)
  where
    a  = poissonMat pp
    b  = poissonVec pp
    triples = cgMethod a b

    -- v:CG法による解, u:厳密解
    (v,_,_) = last triples
    u  = fromListVec [ (pos (i,j), exact_u (pp `ith_x` i) (pp `ith_y` j))
                     | (i,j) <- universe
                     ]
    x  = maxNorm (a <> u .-. b)
    y  = maxNorm (u .-. v)
    n  = fromNat' (sing :: Sing n)

-------------------------------------------------------------------------------
-- Build
-------------------------------------------------------------------------------

buildPoissonDef :: forall n. KnownNat n
                => (Double -> Double -> Double) -- ^ 厳密解
                -> (Double -> Double -> Double) -- ^ Laplacian
                -> PoissonDef n
buildPoissonDef actual_u laplacian = PoissonDef{..}
  where
    rangeX   = (0.0,1.0)
    rangeY   = (0.0,1.0)
    f        = laplacian
    marginX0 = fromListVec [ (i, actual_u (x_ i) y_min ) | i <- universe ]
    margin0Y = fromListVec [ (i, actual_u x_min  (y_ i)) | i <- universe ]
    marginXR = fromListVec [ (i, actual_u (x_ i) y_max ) | i <- universe ]
    marginRY = fromListVec [ (i, actual_u x_max  (y_ i)) | i <- universe ]

    n = fromNat' (sing :: Sing n)
    x_ i = x_min + (x_max-x_min) * (fromIntegral i + 1.0) / (n + 1.0)
    y_ j = y_min + (y_max-y_min) * (fromIntegral j + 1.0) / (n + 1.0)
    (x_min,x_max) = rangeX
    (y_min,y_max) = rangeY


