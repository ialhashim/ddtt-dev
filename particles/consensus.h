// Adapted from code by Samuel Rota Bulò

/**
 * "Robust Region Detection via Consensus Segmentation of Deformable Shapes"
 * E. Rodola, S. Rota Bulo, and D. Cremers
 * Computer Graphics Forum, 33(5), 2014
 *
 * If you use this code, please cite the paper above.
 */
 
#pragma once

#include <unordered_map>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>
using namespace Eigen;

#include "hash_vector.h"

namespace cvlab
{

static inline void hash_keys( const MatrixXd & C, std::vector<double> & in_weights, 
							 MatrixXd & unique_keys, std::vector<double> & out_weights, 
							 std::vector<int> & mapping )
{
	int nt = C.cols();
	int len = C.rows();

	const double* const keys = C.data();
	if(in_weights.size() < nt) in_weights.resize(nt, 1);

	typedef std::unordered_map< std::vector<int>, double > keymap_t;
	keymap_t map;

	std::unordered_map< std::vector<int>, std::vector<int> > key_values;
	mapping.resize(nt, 0);

	std::vector<int> key(len);

	for (int i=0; i<nt; ++i){
		for (int j=0; j<len; ++j)
			key[j] = keys[i*len+j];

		if (map.find(key) == map.end())
			map[key] = in_weights[i];
		else
			map[key] += in_weights[i];

		key_values[key].push_back(i);
	}

	unique_keys = MatrixXd( map.size(), len );
	out_weights.resize( map.size(), 1 );

	int i=0;
	for (keymap_t::const_iterator it=map.begin(); it!=map.end(); ++it, ++i)
	{
		const std::vector<int>& curkey = it->first;
		for (int j=0; j < len; ++j)
			unique_keys( i + (j*map.size()) ) = curkey[j];
		out_weights[i] = it->second;

		// Map back
		{
			auto & cur_key_vals = key_values[curkey];
			for(size_t k = 0; k < cur_key_vals.size(); k++)
				mapping[ cur_key_vals[k] ] = i;
		}
	}
}

static inline void weight_matrix( const std::vector<double> & w, MatrixXd & W )
{
	unsigned int n = (unsigned int)w.size();
	W = MatrixXd::Zero( n, n );
	for(unsigned int i = 0; i < n; i++){
		for(unsigned int j = 0; j < n; j++){
			W(i,j) = w[i] * w[j];
		}
	}
}

static inline void PPC_l2(const MatrixXd & C, const MatrixXd & N, MatrixXd & Y,int k, double toll,int maxIter, int& niter, double& final_err)
{
    int n=C.rows();

	Y = MatrixXd::Zero(k,n);

    /* Y=rand(k,n)+1000; Y=Y./repmat(sum(Y),k,1); */
    Y.setRandom();
    Y.array()+=1000;
    Y.array().rowwise()/=Y.colwise().sum().array(); //column-wise normalization

    /* A=(N.*(Y'Y)-C)+diag(0.5*diag(N)-diag(C)); */
    MatrixXd A=MatrixXd::Zero(n,n);
    A.noalias()=Y.transpose()*Y;
    A=A.array()*N.array()-C.array();
    A.diagonal()=0.5*N.diagonal()-C.diagonal();

    /* G=Y*A */
    MatrixXd G=MatrixXd::Zero(k,n);
    G.noalias()=-Y*A;

    bool notConverged=true;

    niter=0;
    double bestGap;

    while(notConverged && (maxIter<=0 || niter<maxIter)) {
        /* determine search direction */
        int simplexIdx=0;
        int maxIdx=-1, minIdx=-1;
        bestGap=-1;

        for(int i=0; i<n; ++i) {
            int lMaxIdx=0, lMinIdx=-1;
            for(int u=0; u<k; ++u) {
                if(G(u,i)>G(lMaxIdx,i))
                    lMaxIdx=u;
                if(Y(u,i)>0 && (lMinIdx<0 || G(u,i)<G(lMinIdx,i)))
                    lMinIdx=u;
            }
            if(G(lMaxIdx,i)-G(lMinIdx,i)>bestGap) {
                simplexIdx=i;
                maxIdx=lMaxIdx;
                minIdx=lMinIdx;
                bestGap=G(maxIdx,simplexIdx)-G(minIdx,simplexIdx);
            }
        }

        //cout << bestGap << " " << simplexIdx<<" " << maxIdx << " " << minIdx << " " << Y(minIdx,simplexIdx) << endl;

        notConverged=bestGap>toll;
        if(notConverged) {
            /* determine optimal step size */
            double den=(N.row(simplexIdx).array()*((Y.row(maxIdx)-Y.row(minIdx)).array().square())).sum();
            den-=N(simplexIdx,simplexIdx)*(Y(maxIdx,simplexIdx)-Y(minIdx,simplexIdx))*(Y(maxIdx,simplexIdx)-Y(minIdx,simplexIdx));
            den+=N(simplexIdx,simplexIdx)-2*C(simplexIdx,simplexIdx);
            double alpha=den<=0?Y(minIdx,simplexIdx):std::min(bestGap/den,Y(minIdx,simplexIdx));

            //cout << bestGap <<" " << alpha << endl;

            /*update G (undo)*/
            G.noalias()+=Y.col(simplexIdx)*A.row(simplexIdx);

            /*update A*/
            double tmp=A(simplexIdx,simplexIdx);
            A.row(simplexIdx).array()+=alpha*N.row(simplexIdx).array()*(Y.row(maxIdx)-Y.row(minIdx)).array();
            for(int i=0; i<n; ++i)A(i,simplexIdx)=A(simplexIdx,i);
            A(simplexIdx,simplexIdx)=tmp;

            /*update Y*/
            Y(maxIdx,simplexIdx)+=alpha;
            Y(minIdx,simplexIdx)-=alpha;

            /*update G (redo)*/
            G.noalias()-=Y.col(simplexIdx)*A.row(simplexIdx);
            G.col(simplexIdx).noalias()=-Y*A.row(simplexIdx).transpose();

            ++niter;
        }

    }
    //cout << "niter=" << niter << endl;
    final_err = bestGap;
}

// Helper function
MatrixXd reprow( const MatrixXd & X ){
	MatrixXd M = MatrixXd::Zero( X.rows() * X.rows(), X.cols() );
	for(unsigned int i = 0; i < X.rows(); i++){
		size_t idx = (i * X.rows());
		M.block(idx,0,X.rows(),X.cols()) = X.row(i).replicate(X.rows(),1);
	}
	return M;
}
MatrixXd repmat( const MatrixXd & X ){
	MatrixXd M = MatrixXd::Zero( X.rows() * X.rows(), X.cols() );
	for(unsigned int i = 0; i < X.rows(); i++){
		size_t idx = (i * X.rows());
		M.block(idx,0,X.rows(),X.cols()) = X;
	}
	return M;
}

MatrixXd Cbar( const MatrixXd & keys ){
	MatrixXd C = MatrixXd::Zero( keys.rows(), keys.rows() );

	auto X = reprow( keys ), Y = repmat( keys );
	MatrixXd coassoc = MatrixXd::Zero( X.rows(), X.cols() );
	for(unsigned int i = 0; i < X.rows(); i++){
		for(unsigned int j = 0; j < X.cols(); j++){
			if(X(i,j) == Y(i,j))
				coassoc(i,j) = 1;
		}
	}

	VectorXd mean = coassoc.rowwise().mean();
	
	for(unsigned int i = 0; i < keys.rows(); i++){
		for(unsigned int j = 0; j < keys.rows(); j++){
			C(i,j) = (i==j) ? 1.0 : mean( (i * keys.rows()) + j );
		}
	}

	return C;
}

// get the final assignment to clusters as m.a.p.
static inline std::vector<int> assignClusters( const MatrixXd & Y )
{
	std::vector<int> subseg_clusters;

	// Index of maximum element in each a column
	for(unsigned int i = 0; i < Y.cols(); i++){
		VectorXd col = Y.col(i);
		auto max_idx = std::max_element(col.data(), col.data() + Y.rows()) - col.data();
		subseg_clusters.push_back( max_idx );
	}

	return subseg_clusters;
}

struct Consensus{
	double toll, final_err;
	int maxIter, niter;
	std::vector<int> clusters;
};

// Simple interface. Provide a matrix where rows are data points and columns are some assigned clusters.
// Weights are weights for each data point.
static inline Consensus ConsensusClustering(const MatrixXd & clusterings, int k, const MatrixXd & weights = MatrixXd(),
											double toll = 1e-7,int maxIter = 1e5)
{
	int n = clusterings.rows();

	Consensus c;
	c.toll = toll; 
	c.maxIter = maxIter;

	Eigen::MatrixXd cbar = cvlab::Cbar(clusterings);
	Eigen::MatrixXd W = (weights.rows() == n) ? weights : MatrixXd::Ones(n,n);
	Eigen::MatrixXd Y;
	cvlab::PPC_l2(cbar.array() * W.array(), W, Y, k, c.toll, c.maxIter, c.niter, c.final_err);

	c.clusters = assignClusters(Y);

	return c;
}

}
