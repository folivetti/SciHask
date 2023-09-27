{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeApplications #-}
module Main (main) where

import Lib
import Criterion.Main
import System.Random
import Data.Massiv.Array ( computeIO, setComp, S, Comp(..) )
import qualified Numeric.LinearAlgebra as LA
import qualified Data.Massiv.Array as MA
import Data.Massiv.Array.Numeric 

import Numeric.GSL.Minimization (MinimizeMethodD (..), minimizeVD)

g :: StdGen
g = mkStdGen 2020
m :: Int
m = 1000
n :: Int
n = 1000

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
    grad = MA.fromList MA.Seq [grad0, grad1]
    grad0 = MA.sum $ err !*! MA.map dfd0 x
    grad1 = MA.sum $ err !*! MA.map dfd1 x
    dfd0 x0 = -x0 / (t1 + x0)
    dfd1 x0 = t0*x0/(t1 + x0) ^ (2 :: Int)

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
    dfd1 x0 = t0*x0/(t1 + x0) ^ (2 :: Int)
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
    dfd1  = MA.computeAs MA.S . MA.map (\x0 -> t0*x0/(t1 + x0) ^ (2 :: Int))
    dfd00 = MA.computeAs MA.S . MA.map (const 0)
    dfd01 = MA.computeAs MA.S . MA.map (\x0 -> x0 / (t1 + x0) ^ (2 :: Int))
    dfd10 = MA.computeAs MA.S . MA.map (\x0 -> x0/(t1 + x0) ^ (2 :: Int))
    dfd11 = MA.computeAs MA.S . MA.map (\x0 -> -2*t0*x0/(t1 + x0) ^ (3 :: Int))

    t0 = theta MA.! 0
    t1 = theta MA.! 1

thetaPuro :: MassVec
thetaPuro = MA.fromList MA.Seq [200, 0.08]

t0h :: LA.Vector Double 
t0h = LA.fromList [200, 0.08]

setupMult :: IO (MassArray, MassArray, LA.Matrix Double, LA.Matrix Double)
setupMult = do
    let (arr1, arr2) = genTwoArrs m m g
        (arrh1, arrh2) = genTwoH m m g
    pure (arr1, arr2, arrh1, arrh2)

setupDot :: IO (MassArray, MassArray, MassArray, MassArray, MassVec, MassVec, LA.Vector Double, LA.Vector Double)
setupDot = do
    let (arr1, arr2) = genTwoArrs n 1 g
        (arr3, arr4) = genTwoArrs 1 n g
        (v1, v2) = genTwoVecs n g
        (h1, h2) = genTwoHV n g
    pure (arr1, arr2, arr3, arr4, v1, v2, h1, h2)

setupOuter :: IO (MassArray, MassArray, MassVec, MassVec, LA.Vector Double, LA.Vector Double)
setupOuter = do 
    let (arr1, arr2) = genTwoArrs n 1 g 
        (v1, v2) = genTwoVecs n g 
        (h1, h2) = genTwoHV n g 
    pure (arr1, arr2, v1, v2, h1, h2)

setupChol :: IO (MassArray, LA.Herm Double)
setupChol = do 
    let arr = genSymArr m 
        harr = LA.sym (genSymH m)
    pure (arr, harr)

main :: IO ()
main = defaultMain
         [ bgroup "Matrix"
          [env setupMult $ \ ~(arr1, arr2, arrh1, arrh2) ->
             bgroup "main"
               [ bench "Massiv Seq" $ whnfIO (computeIO @S =<< arr1 .><. arr2)
               , bench "Massiv Par" $ whnfIO (computeIO @S =<< arr1 .><. setComp Par arr2)
               , bench "hmatrix" $ whnf (arrh1 LA.<>) arrh2
               ]

          ]
         , bgroup "Dot vecs"
            [env setupDot $ \ ~(arr1, arr2, arr3, arr4, v1, v2, h1, h2) ->
                bgroup "main"
                  [ bench "With transpose" $ whnfIO (computeIO @S =<< MA.computeAs MA.S (MA.transpose arr1) .><. arr2)
                  , bench "With mult transpose" $ whnfIO (computeIO @S =<< MA.multiplyMatricesTransposed arr3 arr4)
                  , bench "!.!" $ whnf (v1 !.!) v2
                  , bench "dotM" $ whnfIO (MA.dotM v1 v2)
                  , bench "dotM Par" $ whnfIO (MA.dotM v1 (setComp Par v2))
                  , bench "hmatrix" $ whnf (LA.dot h1) h2
                  ]
            ]

         , bgroup "Outer vecs"
            [env setupOuter $ \ ~(arr1, arr2, v1, v2, h1, h2) ->
                bgroup "main"
                  [ bench "Massiv Seq" $ whnfIO (computeIO @S =<< arr1 .><. MA.computeAs MA.S (MA.transpose arr2))
                  , bench "Massiv Par" $ whnfIO (computeIO @S =<< arr1 .><. setComp Par (MA.computeAs MA.S (MA.transpose arr2)))
                  , bench "Massiv outer" $ whnfIO (computeIO @S =<< outer v1 (setComp Par v2))
                  , bench "hmatrix" $ whnf (LA.outer h1) h2
                  ]
            ]
         , bgroup "Chol"
            [env setupChol $ \ ~(arr, harr)->
                bgroup "main"
                  [ bench "Massiv Seq" $ whnfIO (computeIO @S =<< cholesky arr)
                  , bench "Massiv Seq" $ whnfIO (computeIO @S =<< unblocked arr)
                  , bench "hmatrix" $ whnf LA.chol harr
                  ]
            ]
         , bgroup "Inv"
            [env setupChol $ \ ~(arr, harr)->
                bgroup "main"
                  [ bench "Massiv Seq" $ whnfIO (computeIO @S =<< invChol arr)
                  --, bench "Massiv Par" $ whnfIO (computeIO @S =<< invChol (setComp Par arr))
                  , bench "hmatrix" $ whnf (\h -> LA.cholSolve (LA.chol h) (LA.ident m)) harr
                  , bench "hmatrix (inv)" $ whnf (LA.inv . LA.unSym) harr
                  ]
            ]
         , bgroup "Opt"
            [ bgroup "main"
                  [ bench "Massiv BFGS" $ whnfIO (computeIO @S =<< minimizeBFGS puro hessPuro 200 1e-16 thetaPuro)
                  , bench "Massiv CG" $ whnfIO (computeIO @S =<< minimizeCG puro 200 1e-16 thetaPuro)
                  , bench "hmatrix BFGS" $ whnf (minimizeVD VectorBFGS2  1e-16 200 1e-1 1e-2 model jacob) t0h
                  , bench "hmatrix CG" $ whnf (minimizeVD ConjugateFR 1e-16 200 1e-1 1e-2 model jacob) t0h

                  , bench "Massiv BFGS - 1 iter" $ whnfIO (computeIO @S =<< minimizeBFGS puro hessPuro 1 1e-16 thetaPuro)
                  , bench "Massiv CG - 1 iter" $ whnfIO (computeIO @S =<< minimizeCG puro 1 1e-16 thetaPuro)
                  , bench "hmatrix BFGS - 1 iter" $ whnf (minimizeVD VectorBFGS2  1e-20 1 1e-1 1e-2 model jacob) t0h
                  , bench "hmatrix CG - 1 iter" $ whnf (minimizeVD ConjugateFR 1e-20 1 1e-1 1e-2 model jacob) t0h
                  ]
            ]
         ]
