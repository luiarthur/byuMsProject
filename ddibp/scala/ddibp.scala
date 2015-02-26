/* K: a list. K[[i]] = the set of dishes owned by customer i.
   l: a vector. l[i] is the cardinality of K[i].
   K.n: the total number of owned dishes. K.n = sum(l).  */

/* C (N x K): a connectivity matrix. C[i,k] = j  =>  customer i connects to cusomter j 
              for dish k.
   c: an ownership vector. c[k] = who owns dish k.  */

/* D (N x N): a distance matrix between all customers. D[i,j] = distance between 
              customers i and j.  D[i,i] = 0. Customer can only connect to previous 
              customers. That is D[1,2] = Inf. But D[2,1] < Inf.  */

/* f: decay function. Maps distance to proximity. High distance => low proximity. 
      f(0) = 1 & f(Inf) = 0.
   A: normalized proximity matrix. A[i,j] = f(D[i,j]) / h[i], where h[i] = sum of 
      the distance b/w i and every other customer.  */



package Rho

// import scala.collection.immutable.Vector.empty
// import breeze.stats.distributions.{NegativeBinomial,Binomial,Gaussian}
// import org.apache.commons.math3.special.Gamma._
import scala.io.Source
import java.io.File
import scala.math.{log,pow,Pi,exp}
import util.Random

class DDIBP {
  def d(Z: Array[Double], a: Double, log: Int): Double = {
    0
  }

  def r(N: Int,a: Double, D: Array[Double]): Array[Double] = {
    Array(0,0)
  }
}
