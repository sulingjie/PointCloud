/*
 * ICP配准算法实现 
 */
#pragma once
#include <iostream>
#include <vector>
#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/time.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>

namespace mypcl {
	// 定义求解最小二乘的方法
	enum {
		QUATERNION,
		FAST, // 该方法已放弃，复现失败
		SVD,
	};
	// 定义计算邻近点的方法
	enum {
		BURTE_FORCE = 0, // this method can not work
		KD_TREE = 1,
	};
	// 定义是否输出调试信息 0:不打印 1:打印
#define SHOW_DEB 0
#define MIN -9999
	using namespace std;

	// print the 3x3 matrix
	void PrintMatrix3(const Eigen::Matrix3d& matrix) {
		printf("    |%6.3f%6.3f%6.3f|\n", matrix(0, 0), matrix(0, 1), matrix(0, 2));
		printf("R = |%6.3f%6.3f%6.3f|\n", matrix(1, 0), matrix(1, 1), matrix(1, 2));
		printf("    |%6.3f%6.3f%6.3f|\n", matrix(2, 0), matrix(2, 1), matrix(2, 2));
	}
	void PrintMatrix4(const Eigen::Matrix4d& matrix) {
		printf("    |%6.3f%6.3f%6.3f%6.3f|\n", matrix(0, 0), matrix(0, 1), matrix(0, 2), matrix(0, 3));
		printf("R = |%6.3f%6.3f%6.3f%6.3f|\n", matrix(1, 0), matrix(1, 1), matrix(1, 2), matrix(1, 3));
		printf("    |%6.3f%6.3f%6.3f%6.3f|\n", matrix(2, 0), matrix(2, 1), matrix(2, 2), matrix(2, 3));
		printf("    |%6.3f%6.3f%6.3f%6.3f|\n", matrix(3, 0), matrix(3, 1), matrix(3, 2), matrix(3, 3));
	}
	// rotation matrix product vector
	vector<float> RotationXVector3(const Eigen::Matrix3f& m, const vector<float>& v) {
		vector<float>vr(3, 0);
		vr[0] = m(0, 0) * v[0] + m(0, 1) * v[1] + m(0, 2) * v[2];
		vr[1] = m(1, 0) * v[0] + m(1, 1) * v[1] + m(1, 2) * v[2];
		vr[2] = m(2, 0) * v[0] + m(2, 1) * v[1] + m(2, 2) * v[2];
		return vr;
	}
	// print vector
	template<typename DataType>
	void PrintVector(const vector<DataType>& A) {
		for (int i = 0; i < A.size(); ++i) {
			cout << A[i] << " ";
		}
		cout << endl;
	}
	template<typename PointSourcePtr, typename PointTargetPtr, typename Scalar = float>
	class ICP {
	private:
		PointSourcePtr cloud_source;
		PointTargetPtr cloud_target;
		int methodMS; // 求解方法 
		int methodCP; // 计算邻近点方法 （穷举法太慢，只能用kd-tree）
		int m_iterations; // 最大迭代次数
		int final_iterations; // 实际迭代次数
		float epsilon; // 两次迭代均方差的差值
		Scalar SquareMeanError; // 均方差
		vector<Scalar> SquareError; // 邻近点之间距离
		Scalar CorrThreshold = 0; // 邻近点阈值，若邻近点之间距离超过该值则忽略
		vector<size_t> Index_corr_src; // 对应点中源点云中的索引
		vector<size_t> Index_corr_tar; // 对应点中目标点云中的索引
		Eigen::Matrix3f rotation_matrix; // rotation matrix
		vector<float> translation_matrix;// translation matrix
		vector<float> unitEV; // quatrenion corrsponding to rotation matrix
	public:
		ICP(const PointSourcePtr source, const PointTargetPtr target) :
			methodMS(1), methodCP(1), unitEV({0,0,0,0}), final_iterations(0){
			cloud_source = source;
			cloud_target = target;
			translation_matrix = vector<float>(3, 0);
		}
		// set the iteration times
		inline void SetIterations(int iterations) { m_iterations = iterations; }
		// set the method sovling Mean Squre errors
		inline void SetSolutionMethod(int method) { methodMS = method; }
		// set the method computing closest points
		inline void SetClosestPointsMethod(int method) { methodCP = method; }
		// set the threshold of distance between correspondences
		template<typename Scalar = float>
		inline void SetCorrThreshold(Scalar threshold) { CorrThreshold = threshold; }
		// set the difference between two iteration
		void SetEpsilon(float e) { epsilon = e; }
		// set the method solving the registration
		void SetMethodMS(int method) { methodMS = method; }
		// set the method finding closest points
		void SetMethodCP(int method) { methodCP = method; }
		// print square mean error
		inline void PrintSME() { cout << "mean square error : " << SquareMeanError << endl; }
		// 向量相乘，两向量均为列向量且项数为3，结果为 A*transpose(B)
		Eigen::Matrix3d PointDot(const pcl::PointXYZ A, const pcl::PointXYZ B) {
			Eigen::Matrix3d matrix = Eigen::Matrix3d::Identity();
			matrix(0, 0) = A.x * B.x;
			matrix(0, 1) = A.x * B.y;
			matrix(0, 2) = A.x * B.z;
			matrix(1, 0) = A.y * B.x;
			matrix(1, 1) = A.y * B.y;
			matrix(1, 2) = A.y * B.z;
			matrix(2, 0) = A.z * B.x;
			matrix(2, 1) = A.z * B.y;
			matrix(2, 2) = A.z * B.z;
			return matrix;
		}
		// vector和数相乘
		template<typename DataType>
		vector<DataType> VectorProduceNum(const vector<DataType>& A, const DataType& num) {
			vector<DataType> B(A.size(), 0);
			for (int i = 0; i < B.size(); ++i) {
				B[i] = A[i] * num;
			}
			return B;
		}
		// 返回迭代次数
		inline int getIterations() { return final_iterations; }
		// 返回最后一次迭代均方差
		inline Scalar getSquareMeanError() { return SquareMeanError; }
		// vector相加
		template <typename DataType>
		vector<DataType> VectorAdd(const vector<DataType>& A, const vector<DataType>& B) {
			if (A.size() != B.size()) {
				cerr << "ERROR:the size of vector is not the same..." << endl;
				return vector<DataType>(0);
			}
			vector<DataType> C(A.size(), 0);
			for (int i = 0; i < A.size(); ++i) {
				C[i] = A[i] + B[i];
			}
			return C;
		}
		// apply the registration to one point
		pcl::PointXYZ Transform(const pcl::PointXYZ& point_src) {
			pcl::PointXYZ point_tar;
			point_tar.x = rotation_matrix(0, 0) * point_src.x + rotation_matrix(0, 1) * point_src.y
				+ rotation_matrix(0, 2) * point_src.z + translation_matrix[0];
			point_tar.y = rotation_matrix(1, 0) * point_src.x + rotation_matrix(1, 1) * point_src.y
				+ rotation_matrix(1, 2) * point_src.z + translation_matrix[1];
			point_tar.z = rotation_matrix(2, 0) * point_src.x + rotation_matrix(2, 1) * point_src.y
				+ rotation_matrix(2, 2) * point_src.z + translation_matrix[2];
			return point_tar;
		}
		// quaternion to rotation matrix
		void QuaterionToRotation(const vector<float>&q, Eigen::Matrix3f &m) {
			float q0 = q[0];float q1 = q[1];float q2 = q[2];float q3 = q[3];
			m(0, 0) = pow(q0, 2) + pow(q1, 2) - pow(q2, 2) - pow(q3, 2);
			m(0, 1) = 2 * (q1 * q2 - q0 * q3);
			m(0, 2) = 2 * (q1 * q3 + q0 * q2);
			m(1, 0) = 2 * (q1 * q2 + q0 * q3);
			m(1, 1) = pow(q0, 2) + pow(q2, 2) - pow(q1, 2) - pow(q3, 2);
			m(1, 2) = 2 * (q2 * q3 - q0 * q1);
			m(2, 0) = 2 * (q1 * q3 - q0 * q2);
			m(2, 1) = 2 * (q2 * q3 + q0 * q1);
			m(2, 2) = pow(q0, 2) + pow(q3, 2) - pow(q1, 2) - pow(q2, 2);
		}
		// 向量相减
		template<typename DataType>
		vector<DataType> VectorSubtract(const vector<DataType>& A, const vector<DataType>& B) {
			if (A.size() != B.size()) {
				cerr << "ERROR:the size of two vector is different..." << endl;
				return vector<DataType>(0);
			}
			vector<DataType>C(A.size(), 0);
			for (int i = 0; i < A.size(); ++i) {
				C[1] = A[i] - B[i];
			}
			return C;
		}
		// get the mean square erroe between source and target
		// 只能在align()后面调用getFitnessScore()
		Scalar getFitnessScore() {
			// find correspondences
			ComputeCorrespondence(cloud_target);
			return SquareMeanError;
		}
		// alian:use four looping steps
		void align(PointTargetPtr& cloud_ICP) {
			*cloud_ICP = *cloud_target;
			int iteration_time = 0;
			if (methodMS == QUATERNION || methodMS == SVD) {
				for (int i = 0; i < m_iterations; ++i) {
					Scalar lastSquareMeanError = SquareMeanError;
					ComputeCorrespondence(cloud_ICP);
#if SHOW_DEB
					PrintSME();
#endif
					if (abs(lastSquareMeanError - SquareMeanError) < epsilon) {
						break;
					}
					++iteration_time;
					ComputeRegistration(cloud_ICP);
					ApplyRegistration(cloud_ICP);
				}
				final_iterations = iteration_time;
#if SHOW_DEB
				cout << "total iteration time : " << iteration_time << endl;
				cout << "the last square mean error : " << SquareMeanError << endl;
#endif
			}
			else if (methodMS == FAST) {
				// 初始化qk,qk-1,qk-2
				vector<float> qk(7, 0), qkm1(7, 0), qkm2(7, 0);
				Scalar dk, dkm1, dkm2;
				// 找到前三次旋转四元数和平移矩阵
				for (int i = 0; i < 3; ++i) {
					ComputeCorrespondence(cloud_ICP);
					ComputeRegistration(cloud_ICP);
					// save the quatarnion and translation matrix
					if (i == 0) {
						qkm2[0] = unitEV[0]; qkm2[1] = unitEV[1]; qkm2[2] = unitEV[2]; qkm2[3] = unitEV[3];
						qkm2[4] = translation_matrix[0]; qkm2[5] = translation_matrix[1]; qkm2[6] = translation_matrix[2];
						dkm2 = SquareMeanError;
						// 
					}
					if (i == 1) {
						qkm1[0] = unitEV[0]; qkm1[1] = unitEV[1]; qkm1[2] = unitEV[2]; qkm1[3] = unitEV[3];
						qkm1[4] = translation_matrix[0]; qkm1[5] = translation_matrix[1]; qkm1[6] = translation_matrix[2];
						dkm1 = SquareMeanError;
					}
					if (i == 2) {
						qk[0] = unitEV[0]; qk[1] = unitEV[1]; qk[2] = unitEV[2]; qk[3] = unitEV[3];
						qk[4] = translation_matrix[0]; qk[5] = translation_matrix[1]; qk[6] = translation_matrix[2];
						dk = SquareMeanError;
					}
					ApplyRegistration(cloud_ICP);
					++iteration_time;
					cout << "the iteration:" << iteration_time << endl;
				}
				cout << "have found the first three num..." << endl;
				float vk, vkm1, vkm2;
				for (int i = 0; i < m_iterations; ++i) {
					// calculate vk,vk-1,vk-2
					cout << "qk : "; PrintVector(qk);
					cout << "qk-1 : "; PrintVector(qkm1);
					cout << "qk-2 : "; PrintVector(qkm2);
					vk = 0;
					vector<float>Dqk = VectorSubtract(qk, qkm1);
					vector<float>Dqkm1 = VectorSubtract(qkm1, qkm2);
					cout << "delta-qk : "; PrintVector(Dqk);
					cout << "delta-qk-1 : "; PrintVector(Dqkm1);

					vkm1 = -sqrt(pow(Dqk[0], 2) + pow(Dqk[1], 2) + pow(Dqk[2], 2) +
						pow(Dqk[3], 2) + pow(Dqk[4], 2) + pow(Dqk[5], 2) + pow(Dqk[6], 2));
					vkm2 = -sqrt(pow(Dqkm1[0], 2) + pow(Dqkm1[1], 2) + pow(Dqkm1[2], 2) +
						pow(Dqkm1[3], 2) + pow(Dqkm1[4], 2) + pow(Dqkm1[5], 2) + pow(Dqkm1[6], 2));
					cout << "vk-1 : " << vkm1 << endl;
					cout << "vk-2 : " << vkm2 << endl;

					// 使用 k 和 k-2 进行线性插值
					float b1 = dk;
					float a1 = (dkm2 - dk) / vkm2;
					// 使用三点进行抛物线插值 checked 
					float a2 = (dkm1 - dk) / pow(vkm1, 2) - ((dkm2 - dk) * pow(vkm1, 2) - (dkm1 - dk) * pow(vkm2, 2)) 
						/ (pow(vkm1, 3) * vkm2 - pow(vkm1, 2) * pow(vkm2, 2));
					float b2 = ((dkm2 - dk) * pow(vkm1, 2) - (dkm1 - dk) * pow(vkm2, 2)) 
						/ (pow(vkm1, 2) * vkm2 - vkm1 * pow(vkm2, 2));
					float c2 = dk;
					// calculate v1 and v2
					float v1 = -b1 / a1;
					float v2 = -b2 / (2 * a2);
					// 计算新的旋转矩阵
					float Dqk_norm = -vkm1;
					float vmax = 25 * Dqk_norm;
					// calculate the next qk
					vector<float> qk_next(7, 0);
					if ((v2 > 0 && v1 > v2 && vmax > v1) || (v2 > 0 && vmax > v2 && v1 > vmax)) {
						qk_next = VectorAdd(qk, VectorProduceNum(Dqk, float(v2 / Dqk_norm)));
					}
					else if ((v1 > 0 && v2 > v1 && vmax > 2) || (v1 > 0 && vmax > v1 && v2 > vmax) || 
						(v2 < 0) || (v1 > 0 && vmax > v1)) {
						qk_next = VectorAdd(qk, VectorProduceNum(Dqk, float(v1 / Dqk_norm)));
					}
					else if ((v1 > vmax) || (v2 > vmax)) {
						qk_next = VectorAdd(qk, VectorProduceNum(Dqk, float(vmax / Dqk_norm)));
					}
					// calculete new rotation matrix and translation matrix
					unitEV[0] = qk_next[0]; unitEV[1] = qk_next[1]; unitEV[2] = qk_next[2]; unitEV[3] = qk_next[3];
					translation_matrix[0] = qk_next[4]; translation_matrix[1] = qk_next[5]; translation_matrix[2] = qk_next[6];
					QuaterionToRotation(unitEV, rotation_matrix);
					// 检查四元数和平移矩阵
					cout << "qk_next : "; PrintVector(qk_next);
					
					// apply rotation and translation matrix
					ApplyRegistration(cloud_ICP);
					// 更新qk
					qkm2 = qkm1;
					qkm1 = qk;
					qk = qk_next;
					++iteration_time;
					cout << "the " << iteration_time << "th iteration" << endl;
				//	return;
				}
			}
		}
		// first step : find the clostest points
		void ComputeCorrespondence(PointTargetPtr& cloud_ICP){
			if (methodCP == BURTE_FORCE) { // 使用穷举法计算邻近点
				// clear the correspondence vector generated at last iteration
				//vector<int>().swap(Index_corr_src);
				//vector<int>().swap(Index_corr_tar);
				//vector<int>().swap(SquareError);
				for (size_t index_source = 0; index_source < cloud_source->points.size(); ++index_source) {
					Scalar minDistance = CorrThreshold;
					size_t minTargetIndex = 0;
					for (size_t index_target = 0; index_target < cloud_target->points.size(); ++index_target) {
						cout << "calculating now... index_source = " << index_source << " index_target = " << index_target << endl;
						//Scalar distance = ComputeDistance3d(
						//	cloud_source->points[index_source].x,
						//	cloud_source->points[index_source].y,
						//	cloud_source->points[index_source].z,
						//	cloud_ICP->points[index_target].x,
						//	cloud_ICP->points[index_target].y,
						//	cloud_ICP->points[index_target].z);
						Scalar test = pow(cloud_source->points[index_source].x - cloud_ICP->points[index_target].x, 2);
						float deltaX = cloud_source->points[index_source].x - cloud_ICP->points[index_target].x;
						float deltaY = cloud_source->points[index_source].y - cloud_ICP->points[index_target].y;
						float deltaZ = cloud_source->points[index_source].z - cloud_ICP->points[index_target].z;
						float distance = pow(float(pow(deltaX, 2) + pow(deltaY, 2) + pow(deltaZ, 2)), 0.5);
						if (distance < minDistance) {
							minDistance = distance;
							minTargetIndex = index_target;
						}
					}
					if (minDistance < CorrThreshold) {
						Index_corr_src.push_back(index_source);
						Index_corr_tar.push_back(minTargetIndex);
						SquareError.push_back(minDistance);
					}
				}
				// calculate the mean square error
				Scalar SquareErrorSum = 0;
				for (size_t i = 0; i < SquareError.size(); ++i) {
					SquareErrorSum += SquareError[i];
				}
				SquareMeanError = SquareErrorSum / SquareError.size();
			}
			else if (methodCP == KD_TREE) { // 使用kdtree 暴力枚举实在是太慢了
				// 由于源点云不变化，所以对源点云建立kdtree
				pcl::KdTreeFLANN<pcl::PointXYZ>kdtree;
				kdtree.setInputCloud(cloud_ICP);
				Index_corr_src.resize(0);
				Index_corr_tar.resize(0);
				SquareError.resize(0);
				for (size_t i = 0; i < cloud_source->points.size(); ++i) {
					pcl::PointXYZ searchPoint;
					searchPoint.x = cloud_source->points[i].x;
					searchPoint.y = cloud_source->points[i].y;
					searchPoint.z = cloud_source->points[i].z;
					int k = 1;
					vector<int> Idx(k);
					vector<Scalar> SD(k);
					if (kdtree.nearestKSearch(searchPoint, k, Idx, SD) > 0) {
						if (SD[0] < CorrThreshold) {
							Index_corr_src.push_back(i);
							Index_corr_tar.push_back(Idx[0]);
							SquareError.push_back(SD[0]);
						}
					}
				/*	if (it == 1) {
						cout << i << " " << Idx[0] << " " << SD[0] << endl;
					}*/
					
				}
				// 计算均方差
				if (SquareError.size() == 0) {
					throw"there si something wrong, SquareError is zero...";
				}
				Scalar SquaerErrorSum = 0;
				for (size_t i = 0; i < SquareError.size(); ++i) {
					SquaerErrorSum += SquareError[i];
				}
				SquareMeanError = SquaerErrorSum / SquareError.size();
			}
		}
		// second step : compute the registration
		void ComputeRegistration(PointTargetPtr& cloud_ICP) {
			if (methodMS == QUATERNION || methodMS == FAST) {
				// 计算源点云和icp点云 邻近点的中心
				pcl::PointXYZ ave_source;
				pcl::PointXYZ ave_ICP;
				CalCentroid(cloud_source, ave_source,Index_corr_src);
				CalCentroid(cloud_ICP, ave_ICP,Index_corr_tar);

#if SHOW_DEB == 1
				cout << "center of source cloud : (" << ave_source.x << "," << ave_source.y
					<< "," << ave_source.z << ")" << endl;
				cout << "center of target cloud : (" << ave_ICP.x << "," << ave_ICP.y 
					<< "," << ave_ICP.z << ")" << endl;
#endif
				// calculate the cross-covariance matrix (3x3)
				Eigen::Matrix3d CroCovarM = Eigen::Matrix3d::Zero();// cross-covariance matrix
				for (size_t i = 0; i < Index_corr_src.size(); ++i) {
					CroCovarM += PointDot(cloud_ICP->points[Index_corr_tar[i]],
						cloud_source->points[Index_corr_src[i]]);
				}
#if SHOW_DEB
				cout << " the temp crocover:" << endl;
				PrintMatrix3(CroCovarM / Index_corr_src.size());
				cout << " the ave_ICP:";
				cout << "(" << ave_ICP.x << "," << ave_ICP.y << "," << ave_ICP.z << ")" << endl;
				cout << "the ave_source:";
				cout << "(" << ave_source.x << "," << ave_source.y << "," << ave_source.z << ")" << endl;
				cout << "the matrix by miu:" << endl;
				PrintMatrix3(PointDot(ave_ICP, ave_source));	
#endif
				CroCovarM = CroCovarM / Index_corr_src.size() - PointDot(ave_ICP, ave_source);
#if SHOW_DEB == 1
				cout << "the size of corr points  : " << Index_corr_src.size() << endl;
				cout << "the cross covariance matrix :" << endl;
				PrintMatrix3(CroCovarM);
#endif
				// calculate the big Q matrix (4x4)
					// first find the vector Delta (3X1)
				Eigen::Matrix3d antiCroCrovarM = CroCovarM - CroCovarM.transpose();
				vector<double> delta(3, 0);
				delta[0] = antiCroCrovarM(1, 2);
				delta[1] = antiCroCrovarM(2, 0);
				delta[2] = antiCroCrovarM(0, 1);
				Eigen::Matrix4d Q = Eigen::Matrix4d::Zero(); // big Q matrix (4x4)
				Q(0, 0) = CroCovarM.trace();
				Q.block<3, 3>(1, 1) = CroCovarM + CroCovarM.transpose() - CroCovarM.trace() * Eigen::Matrix3d::Identity();
				Q(1, 0) = delta[0]; Q(2, 0) = delta[1]; Q(3, 0) = delta[2];
				Q(0, 1) = delta[0]; Q(0, 2) = delta[1]; Q(0, 3) = delta[2];
#if SHOW_DEB == 1
				cout << "the big Q matricx:" << endl;
				PrintMatrix4(Q);
#endif
				// capture the unit eigenvector corresponding to the maximun eigenvalue of matrix Q
				Eigen::EigenSolver<Eigen::Matrix4d>ev(Q,true);
				    // contract the index of the maximum eigenvalue
				float maxEV = MIN;
				int IndMaxEV = -1;
				for (int i = 0; i < 4; ++i) {
					float tempEV = ev.eigenvalues()[i].real();
					if (tempEV > maxEV) {
						maxEV = tempEV;
						IndMaxEV = i;
					}
				}
					// calculate the unit eigenvector
				
				unitEV[0] = ev.eigenvectors()(0,IndMaxEV).real();
				unitEV[1] = ev.eigenvectors()(1,IndMaxEV).real();
				unitEV[2] = ev.eigenvectors()(2,IndMaxEV).real();
				unitEV[3] = ev.eigenvectors()(3,IndMaxEV).real();
				float ModUnitEV = sqrt(pow(unitEV[0], 2) + pow(unitEV[1], 2) + pow(unitEV[2], 2) + pow(unitEV[3], 2));
				unitEV[0] /= ModUnitEV; unitEV[1] /= ModUnitEV; unitEV[2] /= ModUnitEV; unitEV[3] /= ModUnitEV;
#if SHOW_DEB
				cout << "the eigenvalue: " << endl << ev.eigenvalues() << endl;
				cout << "the eigenvector: " << endl << ev.eigenvectors() << endl;
				cout << "the index of max eigenvalue: " << IndMaxEV << endl;
				cout << "the corresponding eigenvector : ";
				for (int i = 0; i < unitEV.size(); ++i) {
					cout << unitEV[i] << " ";
				}
				cout << endl;
#endif
				// calculate the corresponding rotation matrix and translation matrix
				QuaterionToRotation(unitEV, rotation_matrix);
#if SHOW_DEB
				cout << "rotation acquired by paper:" << endl;
				PrintMatrix3(rotation_matrix.cast<double>());
#endif
				vector<float> miuX(3, 0);
				vector<float> miuP(3, 0);
				miuX[0] = ave_source.x; miuX[1] = ave_source.y; miuX[2] = ave_source.z;
				miuP[0] = ave_ICP.x; miuP[1] = ave_ICP.y; miuP[2] = ave_ICP.z;
				vector<float>temp(3, 0);
				temp = RotationXVector3(rotation_matrix, miuP);
				
				translation_matrix[0] = miuX[0] - temp[0];
				translation_matrix[1] = miuX[1] - temp[1];
				translation_matrix[2] = miuX[2] - temp[2];
			}
			else if(methodMS == SVD){
				/********* calculate the centroid of clouds *****/
				pcl::PointXYZ ave_source;
				pcl::PointXYZ ave_ICP;
				Eigen::Vector3d ave_source_vec;
				Eigen::Vector3d ave_ICP_vec;
				CalCentroid(cloud_source, ave_source,Index_corr_src);
				CalCentroid(cloud_ICP, ave_ICP,Index_corr_tar);
				ave_source_vec(0) = ave_source.x; ave_source_vec(1) = ave_source.y; ave_source_vec(2) = ave_source.z;
				ave_ICP_vec(0) = ave_ICP.x; ave_ICP_vec(1) = ave_ICP.y; ave_ICP_vec(2) = ave_ICP.z;
				cout << ave_source_vec << endl;
				cout << ave_ICP_vec << endl;
				// calculate the new clouds
				pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_source_bar(new pcl::PointCloud<pcl::PointXYZ>);
				pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ICP_bar(new pcl::PointCloud<pcl::PointXYZ>);
				cloud_source_bar->resize(cloud_source->size());
				cloud_ICP_bar->resize(cloud_ICP->size());
				// calculate the location 
				for (size_t i = 0; i < cloud_source->size(); ++i) {
					cloud_source_bar->points[i].x = cloud_source->points[i].x - ave_source.x;
					cloud_source_bar->points[i].y = cloud_source->points[i].y - ave_source.y;
					cloud_source_bar->points[i].z = cloud_source->points[i].z - ave_source.z;
				}
				for (size_t i = 0; i < cloud_ICP->size(); ++i) {
					cloud_ICP_bar->points[i].x = cloud_ICP->points[i].x - ave_ICP.x;
					cloud_ICP_bar->points[i].y = cloud_ICP->points[i].y - ave_ICP.y;
					cloud_ICP_bar->points[i].z = cloud_ICP->points[i].z - ave_ICP.z;
				}
				/********* calculate the 3x3 matrix H*/
				Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
				for (size_t i = 0; i < Index_corr_src.size(); ++i) {
					Eigen::Vector3d p_source;
					p_source(0) = cloud_source_bar->points[Index_corr_src[i]].x;
					p_source(1) = cloud_source_bar->points[Index_corr_src[i]].y;
					p_source(2) = cloud_source_bar->points[Index_corr_src[i]].z;
					Eigen::Vector3d p_ICP;
					p_ICP(0) = cloud_ICP_bar->points[Index_corr_tar[i]].x;
					p_ICP(1) = cloud_ICP_bar->points[Index_corr_tar[i]].y;
					p_ICP(2) = cloud_ICP_bar->points[Index_corr_tar[i]].z;

					// flush matrix H
					H = H + p_ICP * p_source.transpose();
				}
				// SVD of H
				Eigen::JacobiSVD<Eigen::Matrix3d> H_svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
				// calculate the 3*3 matrix X
				Eigen::Matrix3d X = H_svd.matrixV() * H_svd.matrixU().transpose();
				//cout << X << endl;
				//cout << "det of X: " << X.determinant() << endl;
				if (X.determinant() > 0) {
					Eigen::Matrix3d rotation = X;
					Eigen::Vector3d translation = ave_source_vec - rotation * ave_ICP_vec;
					cout << translation << endl;
					rotation_matrix = rotation.cast<float>();
					translation_matrix[0] = translation(0);
					translation_matrix[1] = translation(1);
					translation_matrix[2] = translation(2);
				}
				else {
					cout << "failed, matrix X is positive..." << endl;
					abort();
				}

			}
		}
		// thrid step : apply the registration to aquare new cloud
		void ApplyRegistration(PointTargetPtr& cloud_ICP) {
			for (size_t i = 0; i < cloud_ICP->points.size(); ++i) {
				cloud_ICP->points[i] = Transform(cloud_ICP->points[i]);
			}
		}
		/******************** middle function***************/
		// calculate the controid of a set of point cloud
		template<typename PointCloudType,typename PointType>
		void CalCentroid(const PointCloudType& cloud, PointType& ave_point,const vector<size_t>& Index_corr) {
			PointType ave_source;
			Scalar tempX = 0;
			Scalar tempY = 0;
			Scalar tempZ = 0;
			for (size_t i = 0; i < Index_corr.size(); ++i) {
				tempX += cloud->points[Index_corr[i]].x;
				tempY += cloud->points[Index_corr[i]].y;
				tempZ += cloud->points[Index_corr[i]].z;
			}
			ave_source.x = tempX / Index_corr.size();
			ave_source.y = tempY / Index_corr.size();
			ave_source.z = tempZ / Index_corr.size();
			ave_point = ave_source;
		}
	};
};


