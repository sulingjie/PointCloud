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
		int _cloudsNum; //��������
		int _iterNum;// ����������
		float _corrThreshold;//�ڽ�����ֵ
		float _alpha;// ��������
		float _meanSquaredError;// ����������ƽ������
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
		void FindCorrespondence();// ʹ�� _clouds_in �� _clouds_ICP �����ڽ���
		void UpdateCloudsPosition(); // ��R��TӦ�õ� _clouds_ICP 
	public:
		GloballyICP(const vector<pcl::PointCloud<pcl::PointXYZ>>& point_clouds, int Ifix);
		void align(vector<pcl::PointCloud<pcl::PointXYZ>>& clouds_out);// ������׼����
		void setIterationNum(int iterNum) { _iterNum = iterNum; }
		void setCorrThreshold(float threshold) { _corrThreshold = threshold; }
		void setAlpha(float alpha) { _alpha = alpha; }
	};
}