{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeApplications #-}
import Lib
import Data.Massiv.Array ( computeIO, Ix2(..), S )
import qualified Numeric.LinearAlgebra as LA
import qualified Data.Massiv.Array as MA
import Data.Massiv.Array.Numeric
import Numeric.GSL.Minimization (MinimizeMethodD (..), minimizeVD)

m :: Int
m = 10

nonlinEx :: MassVec -> (Double, MassVec)
nonlinEx theta = (f^(2 :: Int), MA.fromList MA.Seq [1*2*f, 4 * t1 * 2 * f, 9*2*f*t2^(2 :: Int)])
  where
    f  = t0 + 2 * t1^(2 :: Int) + 3 * t2^(3 :: Int)
    t0 = theta MA.! 0
    t1 = theta MA.! 1
    t2 = theta MA.! 2

hessEx :: MassVec -> MassArray
hessEx theta = MA.makeArray MA.Seq (MA.Sz (3 :. 3)) f
  where
    g  = t0 + 2 * t1^(2 :: Int) + 3 * t2^(3 :: Int)

    f :: MA.Ix2 -> Double
    f (0 :. 0) = 2
    f (0 :. 1) = 8 * t1
    f (0 :. 2) = 18 * t2^(2 :: Int)
    f (1 :. 0) = 8 * t1
    f (1 :. 1) = 8*g + 8*t1*4*t1
    f (1 :. 2) = 8*t1*9*t2^(2 :: Int)
    f (2 :. 0) = 18 * t2^(2 :: Int)
    f (2 :. 1) = 18*t2^(2 :: Int)*4*t1
    f (2 :. 2) = 36*g*t2 + 18*t2^(2 :: Int)*9*t2^(2 :: Int)
    f _        = 0

    t0 = theta MA.! 0
    t1 = theta MA.! 1
    t2 = theta MA.! 2

thetaEx :: MassVec
thetaEx = MA.fromList MA.Seq [1, 2, 3]

puro :: MassVec -> (Double, MassVec)
puro theta = (gauss, grad)
  where
    x, y :: MassVec 
    x = MA.fromList MA.Seq [0.02, 0.02, 0.06, 0.06, 0.11, 0.11, 0.22, 0.22, 0.56, 0.56, 1.10, 1.10]
    y = MA.fromList MA.Seq [76, 47, 97, 107, 123, 139, 159, 152, 191, 201, 207, 200]
    ypred = MA.map f x
    err = MA.zipWith (-) y ypred
    gauss = (0.5 *) . MA.sum $ err !*! err
    f x0 = t0*x0 / (t1 + x0)
    -- grad = 0.5 sum (err^(2 :: Int)) = sum err * err' = sum err * (-ypred')
    -- hess = sum (ypred'^(2 :: Int) - err * ypred'')
    grad = MA.fromList MA.Seq [grad0, grad1]
    grad0 = MA.sum $ err !*! MA.map dfd0 x
    grad1 = MA.sum $ err !*! MA.map dfd1 x
    dfd0 x0 = -x0 / (t1 + x0)
    dfd1 x0 = t0*x0/(t1 + x0)^(2 :: Int)

    t0 = theta MA.! 0
    t1 = theta MA.! 1

model :: LA.Vector Double -> Double 
model theta = gauss
  where
    x = LA.fromList [0.02, 0.02, 0.06, 0.06, 0.11, 0.11, 0.22, 0.22, 0.56, 0.56, 1.10, 1.10]
    y = LA.fromList [76, 47, 97, 107, 123, 139, 159, 152, 191, 201, 207, 200]

    ypred = LA.cmap f x
    err = y - ypred
    gauss = (0.5 *) . LA.sumElements $ err ^ (2 :: Int)
    f x0 = t0*x0 / (t1 + x0)

    t0 = theta LA.! 0
    t1 = theta LA.! 1

jacob :: LA.Vector Double -> LA.Vector Double 
jacob theta = LA.fromList [grad0, grad1] 
  where
    x = LA.fromList [0.02, 0.02, 0.06, 0.06, 0.11, 0.11, 0.22, 0.22, 0.56, 0.56, 1.10, 1.10]
    y = LA.fromList [76, 47, 97, 107, 123, 139, 159, 152, 191, 201, 207, 200]

    ypred = LA.cmap f x
    err = y - ypred
    f x0 = t0*x0 / (t1 + x0)

    grad0 = LA.sumElements $ err * LA.cmap dfd0 x
    grad1 = LA.sumElements $ err * LA.cmap dfd1 x
    dfd0 x0 = -x0 / (t1 + x0)
    dfd1 x0 = t0*x0/(t1 + x0)^(2 :: Int)
    t0 = theta LA.! 0
    t1 = theta LA.! 1

hessPuro :: MassVec -> MassArray
hessPuro theta = MA.fromLists' MA.Seq [[hess00, hess01], [hess10, hess11]]
  where
    x, y :: MassVec
    x = MA.fromList MA.Seq [0.02, 0.02, 0.06, 0.06, 0.11, 0.11, 0.22, 0.22, 0.56, 0.56, 1.10, 1.10]
    y = MA.fromList MA.Seq [76, 47, 97, 107, 123, 139, 159, 152, 191, 201, 207, 200]
    ypred = MA.map f x
    err :: MassVec
    err = MA.computeAs MA.S $ MA.zipWith (-) y ypred
    f x0 = t0*x0 / (t1 + x0)

    hess00, hess01, hess10, hess11 :: Double
    hess00 = MA.sum $ (dfd0 x !*! dfd0 x) !-! (err !*! dfd00 x)
    hess01 = MA.sum $ (dfd0 x !*! dfd1 x) !-! (err !*! dfd01 x)
    hess10 = MA.sum $ (dfd1 x !*! dfd0 x) !-! (err !*! dfd10 x)
    hess11 = MA.sum $ (dfd1 x !*! dfd1 x) !-! (err !*! dfd11 x)

    dfd0, dfd1, dfd00, dfd01, dfd10, dfd11 :: MassVec -> MassVec
    dfd0  = MA.computeAs MA.S . MA.map (\x0 -> -x0 / (t1 + x0))
    dfd1  = MA.computeAs MA.S . MA.map (\x0 -> t0*x0/(t1 + x0)^(2 :: Int))
    dfd00 = MA.computeAs MA.S . MA.map (const 0)
    dfd01 = MA.computeAs MA.S . MA.map (\x0 -> x0 / (t1 + x0)^(2 :: Int))
    dfd10 = MA.computeAs MA.S . MA.map (\x0 -> x0/(t1 + x0)^(2 :: Int))
    dfd11 = MA.computeAs MA.S . MA.map (\x0 -> -2*t0*x0/(t1 + x0)^(3 :: Int))

    t0 = theta MA.! 0
    t1 = theta MA.! 1

thetaPuro :: MassVec
thetaPuro = MA.fromList MA.Seq [200, 0.08]

t0h :: LA.Vector Double 
t0h = LA.fromList [200, 0.08]

isApproxEqualMM :: MassArray -> MassArray -> Bool
isApproxEqualMM arr1 arr2 = all (\(x,y) -> abs (x-y) < 1e-6) $ zip (MA.toList arr1) (MA.toList arr2)
isApproxEqualMH :: MassArray -> LA.Matrix Double -> Bool
isApproxEqualMH arr h = all (\(x,y) -> abs (x-y) < 1e-6) $ zip (MA.toList arr) (LA.toList $ LA.flatten h)

main :: IO ()
main = do
  let arr   = genSymArr m
      arr2  = genSymArr m
      h     = LA.sym $ genSymH m
      cholH = LA.tr' $ LA.chol h
      invH  = LA.cholSolve (LA.chol h) (LA.ident m)
  cholM  <- computeIO @S =<< cholesky arr
  cholM' <- computeIO @S =<< cholesky arr
  invM   <- invChol arr



  if (isApproxEqualMH cholM cholH)
     then putStrLn "cholesky - hmatrix...ok"
     else do putStrLn "cholesky - hmatrix failed:"
             print cholM
             print cholH
  if (isApproxEqualMH cholM' cholH)
     then putStrLn "unblocked - hmatrix...ok"
     else do putStrLn "unblocked - hmatrix failed:"
             print cholM'
             print cholH
  if (isApproxEqualMH cholM cholH)
     then putStrLn "invChol - hmatrix...ok"
     else do putStrLn "invChol - hmatrix failed:"
             print invM
             print invH

  if (isApproxEqualMM arr arr2)
     then putStrLn "SymArr matches...ok"
     else do putStrLn "SymArr test failed:"
             print arr
             print arr2
  minimizeBFGS nonlinEx hessEx 100 1e-16 thetaEx >>= print
  minimizeCG nonlinEx 100 1e-16 thetaEx >>= print
  putStr "BFGS " >> minimizeBFGS puro hessPuro 200 1e-16 thetaPuro >>= print
  print thetaPuro >> putStr "CG " >> minimizeCG puro 200 1e-16 thetaPuro >>= print
  let  (opt1, _) = minimizeVD ConjugateFR 1e-16 200 1e-1 1e-2 model jacob t0h
  putStr "hmatrix CG " >> print opt1
  let  (opt2, _) = minimizeVD VectorBFGS2 1e-16 200 1e-1 1e-2 model jacob t0h
  putStr "hmatrix BFGS2 " >> print opt2
