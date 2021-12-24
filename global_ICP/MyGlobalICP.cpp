#include "MyGlobalICP.h"
#include <iostream>
#include <vector>
#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/time.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>

using namespace std;
#define DEBUG 0
namespace mypcl {
	GloballyICP::GloballyICP(const vector<pcl::PointCloud<pcl::PointXYZ>>& point_clouds, int Ifix) {
		if (point_clouds.size() == 0 || Ifix >= point_clouds.size()) {
			throw"size of point clouds is 0 or Ifix is too large";
		}
		_cloudsNum = point_clouds.size();
		_clouds_in = point_clouds;
		_clouds_ICP = _clouds_in;
		_Qclouds_ICP.resize(_cloudsNum);
		for (int i = 0; i < _cloudsNum; ++i) {
			_Qclouds_ICP[i].resize(_clouds_in[i].size());
		}
		Eigen::Matrix3f tempRotationMatrix = Eigen::Matrix3f::Identity();
		Eigen::Vector3f tempTransformation = Eigen::Vector3f::Zero();
		pcl::PointXYZ tempCentroid(0, 0, 0);
		_correspondence_index.resize(_cloudsNum);
		_correspondenceNum.resize(_cloudsNum);
		_meanError.resize(_cloudsNum);
		for (int i = 0; i < _cloudsNum; ++i) {
			_Rotation_Matrices.push_back(tempRotationMatrix);
			_delta_Rotation_Matrices.push_back(tempRotationMatrix);
			_Transformaton_Matrices.push_back(tempTransformation);
			_Centroids.push_back(tempCentroid);
			vector<size_t>tempCorr(0);
			_correspondence_index[i].resize(_cloudsNum);
			_correspondenceNum[i].resize(_cloudsNum);
			_meanError[i].resize(_cloudsNum);
			for (int j = 0; j < _cloudsNum; ++j) {
				_correspondence_index[i][j].resize(0);
			}
		}
		_Ifix = Ifix;
		_iterNum = 50; // 默认最大迭代次数为50
		_corrThreshold = 3.0; // 默认邻近点阈值为3
		_alpha = 0.0001;// 迭代步长默认为1
	}
	void GloballyICP::align(vector<pcl::PointCloud<pcl::PointXYZ>>& clouds_out) {
		for (int indexIter = 0; indexIter < _iterNum; ++indexIter) {
			cout << "iteration : " << indexIter << endl;
			FindCentroid();
			for (size_t indexCloud = 0; indexCloud < _cloudsNum; ++indexCloud) {
				int centroidX = _Centroids[indexCloud].x;
				int centroidY = _Centroids[indexCloud].y;
				int centroidZ = _Centroids[indexCloud].z;
				for (size_t indexPoint = 0; indexPoint < _clouds_ICP[indexCloud].points.size(); ++indexPoint) {
					_Qclouds_ICP[indexCloud][indexPoint].x = _clouds_ICP[indexCloud][indexPoint].x - centroidX;
					_Qclouds_ICP[indexCloud][indexPoint].y = _clouds_ICP[indexCloud][indexPoint].y - centroidY;
					_Qclouds_ICP[indexCloud][indexPoint].z = _clouds_ICP[indexCloud][indexPoint].z - centroidZ;
				}
			}
			/* step one : find correspondences */
			FindCorrespondence();
			cout << "MSE : " << _meanSquaredError << endl;
			/* step two : 梯度法求解旋转矩阵R */
			// 遍历计算梯度
			for (int indexR = 0; indexR < _cloudsNum; ++indexR) {
				//cout << indexR << endl;
				if (indexR == _Ifix) {
					_delta_Rotation_Matrices[indexR] = Eigen::Matrix3f::Zero();
					continue;
				}
				Eigen::Matrix3f tempDeltaRotation = Eigen::Matrix3f::Zero();
				int pointNum = 0;
				for (int indexCloud = 0; indexCloud < _cloudsNum; ++indexCloud) {
					if (indexCloud == indexR)
						continue;
					vector<size_t>correspondenceM = _correspondence_index[indexR][indexCloud];
					vector<size_t>correspondenceN = _correspondence_index[indexCloud][indexR];
					for (size_t indexPoint = 0; indexPoint < _correspondenceNum[indexR][indexCloud]; ++indexPoint) {
						//cout << _Qclouds_ICP[indexR].points[correspondenceM[indexPoint]] << endl;
						Eigen::Vector3f Mpoint(_Qclouds_ICP[indexR].points[correspondenceM[indexPoint]].x,
							_Qclouds_ICP[indexR].points[correspondenceM[indexPoint]].y,
							_Qclouds_ICP[indexR].points[correspondenceM[indexPoint]].z);
						//cout << _Qclouds_ICP[indexCloud].points[correspondenceN[indexPoint]] << endl;
						Eigen::Vector3f Npoint(_Qclouds_ICP[indexCloud].points[correspondenceN[indexPoint]].x,
							_Qclouds_ICP[indexCloud].points[correspondenceN[indexPoint]].y,
							_Qclouds_ICP[indexCloud].points[correspondenceN[indexPoint]].z);
						//cout << "Mpoint:" << endl << Mpoint << endl;
						//cout << "Npoint:" << endl << Npoint << endl;
						tempDeltaRotation += 2 * (_Rotation_Matrices[indexR] * Mpoint * Mpoint.transpose() -
							_Rotation_Matrices[indexCloud] * Npoint * Mpoint.transpose());
						++pointNum;
					}
				}
				if (pointNum == 0) {
					_delta_Rotation_Matrices[indexR] = Eigen::Matrix3f::Zero();
				}
				else {
					_delta_Rotation_Matrices[indexR] = tempDeltaRotation / pointNum;
				}
				//cout << "delta rotation Matrix:" << indexR << endl << _delta_Rotation_Matrices[indexR] << endl;
			}
			// 计算更新后的旋转矩阵 R
			for (int indexR = 0; indexR < _cloudsNum; ++indexR) {
				_Rotation_Matrices[indexR] = _Rotation_Matrices[indexR] - _alpha * _delta_Rotation_Matrices[indexR];
			}
			//cout << "RotationMatrix:" << endl << _Rotation_Matrices[1] << endl;
			// 计算平移矩阵
			Eigen::Vector3f centroidF(_Centroids[_Ifix].x, _Centroids[_Ifix].y, _Centroids[_Ifix].z);
			for (int indexT = 0; indexT < _cloudsNum; ++indexT) {
				if (indexT == _Ifix)
					_Transformaton_Matrices[indexT] = Eigen::Vector3f::Zero();
				else {
					Eigen::Vector3f centroidM(_Centroids[indexT].x, _Centroids[indexT].y, _Centroids[indexT].z);
					_Transformaton_Matrices[indexT] = centroidF - _Rotation_Matrices[indexT] * centroidM;
					//cout << "Transformation matrix " << indexT << " : " << endl;
					//cout << _Transformaton_Matrices[indexT] << endl;
				}
			}
			/* step3 : 将R，T应用到 _clouds_ICP  */
			UpdateCloudsPosition();
		}
		clouds_out = _clouds_ICP;
	}
	void GloballyICP::FindCentroid() {
		for (size_t indexCloud = 0; indexCloud < _cloudsNum; ++indexCloud) {
			pcl::PointXYZ centroid(0, 0, 0);
			for (size_t indexPoint = 0; indexPoint < _clouds_ICP[indexCloud].points.size(); ++indexPoint) {
				centroid.x += _clouds_ICP[indexCloud].points[indexPoint].x;
				centroid.y += _clouds_ICP[indexCloud].points[indexPoint].y;
				centroid.z += _clouds_ICP[indexCloud].points[indexPoint].z;
			}
			centroid.x = centroid.x / _clouds_ICP[indexCloud].points.size();
			centroid.y = centroid.y / _clouds_ICP[indexCloud].points.size();
			centroid.z = centroid.z / _clouds_ICP[indexCloud].points.size();
			_Centroids[indexCloud] = centroid;
		}
#if DEBUG == 1
		cout << "centroids:" << endl;
		for (size_t i = 0; i < _cloudsNum; ++i)
			cout << _Centroids[i] << endl;
#endif
	}
	void GloballyICP::FindCorrespondence() {
		// 清空 _correspondence_index, _correspondenceNum, _meanError
		for (int i = 0; i < _cloudsNum; ++i) {
			for (int j = 0; j < _cloudsNum; ++j) {
				_correspondenceNum[i][j] = 0;
				_meanError[i][j] = 0;
				_correspondence_index[i][j].clear();
			}
		}
		size_t totalNum = 0;
		_meanSquaredError = 0;
		for (int i = 0; i < _cloudsNum; ++i) {
			vector<size_t>a(0);
			_correspondence_index[i][i] = a;
			_correspondenceNum[i][i] = 0;
			_meanError[i][i] = 0;
			pcl::KdTreeFLANN<pcl::PointXYZ>kdtree;
			kdtree.setInputCloud(_clouds_ICP[i].makeShared());
			for (int j = i + 1; j < _cloudsNum; ++j) {
				double isumError = 0;
				for (size_t indexPoint = 0; indexPoint < _clouds_ICP[j].size(); ++indexPoint) {
					pcl::PointXYZ searchPoint;
					searchPoint.x = _clouds_ICP[j].points[indexPoint].x;
					searchPoint.y = _clouds_ICP[j].points[indexPoint].y;
					searchPoint.z = _clouds_ICP[j].points[indexPoint].z;
					int k = 1;
					vector<int> Idx(k);
					vector<float> SD(k);
					if (kdtree.nearestKSearch(searchPoint, k, Idx, SD) > 0) {
						if (SD[0] < _corrThreshold*_corrThreshold) {
							_correspondence_index[j][i].push_back(indexPoint);
							_correspondence_index[i][j].push_back(Idx[0]);
							isumError += SD[0];
							_meanSquaredError += sqrt(SD[0]);
							++totalNum;
						}
					}
				}
				_correspondenceNum[i][j] = _correspondence_index[i][j].size();
				_correspondenceNum[j][i] = _correspondence_index[j][i].size();
				if (_correspondenceNum[i][j] == 0)
					_meanError[i][j] = _meanError[j][i] = 0;
				else {
					_meanError[i][j] = isumError / _correspondenceNum[i][j];
					_meanError[j][i] = isumError / _correspondenceNum[j][i];
				}				
			}
		}
		_meanSquaredError = _meanSquaredError / totalNum;
#if DEBUG == 1
		cout << "correspondence num:" << endl;
		for (int i = 0; i < _correspondenceNum.size(); ++i) {
			for (int j = 0; j < _correspondenceNum.size(); ++j) {
				cout << _correspondenceNum[i][j] << " ";
			}
			cout << endl;
		}
		cout << "sum error:" << endl;
		for (int i = 0; i < _meanError.size(); ++i) {
			for (int j = 0; j < _meanError.size(); ++j) {
				cout << _meanError[i][j] << " ";
			}
			cout << endl;
		}
#endif
	}
	void GloballyICP::UpdateCloudsPosition() {
		Eigen::Vector3f tempPoint = Eigen::Vector3f::Ones();
		for (int indexCloud = 0; indexCloud < _cloudsNum; ++indexCloud) {
			for (size_t indexPoint = 0; indexPoint < _clouds_ICP[indexCloud].size(); ++indexPoint) {
				//cout << "point before transform:" << endl;
				//cout << _clouds_ICP[indexCloud][indexPoint] << endl;
				tempPoint(0) = _clouds_ICP[indexCloud][indexPoint].x;
				tempPoint(1) = _clouds_ICP[indexCloud][indexPoint].y;
				tempPoint(2) = _clouds_ICP[indexCloud][indexPoint].z;
				tempPoint = _Rotation_Matrices[indexCloud] * tempPoint + _Transformaton_Matrices[indexCloud];
				//cout << "tempPoint:" << endl << tempPoint << endl;
				_clouds_ICP[indexCloud][indexPoint].x = tempPoint(0);
				_clouds_ICP[indexCloud][indexPoint].y = tempPoint(1);
				_clouds_ICP[indexCloud][indexPoint].z = tempPoint(2);
	/*			cout << "point after transform:" << endl;
				cout << _clouds_ICP[indexCloud][indexPoint] << endl;*/
			}
			/*for (size_t indexPoint = 0; indexPoint < _Qclouds_ICP[indexCloud].size(); ++indexPoint) {
				tempPoint(0) = _Qclouds_ICP[indexCloud][indexPoint].x;
				tempPoint(1) = _Qclouds_ICP[indexCloud][indexPoint].y;
				tempPoint(2) = _Qclouds_ICP[indexCloud][indexPoint].z;
				tempPoint = _Rotation_Matrices[indexCloud] * tempPoint;
				_Qclouds_ICP[indexCloud][indexPoint].x = tempPoint(0);
				_Qclouds_ICP[indexCloud][indexPoint].y = tempPoint(1);
				_Qclouds_ICP[indexCloud][indexPoint].z = tempPoint(2);
			}*/
		}
	}
}