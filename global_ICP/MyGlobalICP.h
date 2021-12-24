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

using namespace std;
namespace mypcl {
	
	class GloballyICP {
	private:
		int _cloudsNum; //点云数量
		int _iterNum;// 最大迭代次数
		float _corrThreshold;//邻近点阈值
		float _alpha;// 迭代步长
		float _meanSquaredError;// 相邻两点间的平均距离
		vector<pcl::PointCloud<pcl::PointXYZ>> _clouds_in;
		vector<pcl::PointCloud<pcl::PointXYZ>> _clouds_ICP;
		vector<pcl::PointCloud<pcl::PointXYZ>> _clouds_out;
		vector<pcl::PointCloud<pcl::PointXYZ>> _Qclouds_ICP;
		vector<vector<vector<size_t>>>_correspondence_index;
		vector<vector<size_t>>_correspondenceNum;
		vector<vector<double>>_meanError;
		vector<Eigen::Matrix3f>_Rotation_Matrices;
		vector<Eigen::Matrix3f>_delta_Rotation_Matrices;
		vector<Eigen::Vector3f>_Transformaton_Matrices;
		vector<pcl::PointXYZ>_Centroids; 
		int _Ifix;
		void FindCentroid();
		void FindCorrespondence();// 使用 _clouds_in 和 _clouds_ICP 计算邻近点
		void UpdateCloudsPosition(); // 将R，T应用到 _clouds_ICP 
	public:
		GloballyICP(const vector<pcl::PointCloud<pcl::PointXYZ>>& point_clouds, int Ifix);
		void align(vector<pcl::PointCloud<pcl::PointXYZ>>& clouds_out);// 核心配准函数
		void setIterationNum(int iterNum) { _iterNum = iterNum; }
		void setCorrThreshold(float threshold) { _corrThreshold = threshold; }
		void setAlpha(float alpha) { _alpha = alpha; }
	};
}