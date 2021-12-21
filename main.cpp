#include <iostream>
#include <vector>
#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/time.h>
#include "MyICP.h"

#define FILTER 0
using namespace std;
int main(int argc, char** argv) {
	//Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
	//Eigen::Vector3d a;
	//a << 1,2, 1;
	//cout << a * a.transpose() << endl;
	//M = a * a.transpose();
	//cout << M<<endl;
	//cout << a * a.transpose() + a * a.transpose() << endl;
	//Eigen::JacobiSVD<Eigen::Matrix3d> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
	//cout << svd.matrixU() << endl;
	//cout << svd.matrixV() << endl;
	//cout << svd.singularValues() <<endl;
	//return 0;
#if FILTER
	pcl::PCLPointCloud2::Ptr cloud (new pcl::PCLPointCloud2);
	pcl::PCLPointCloud2::Ptr cloud_filtered (new pcl::PCLPointCloud2);
		// 填入点云数据
	pcl::PLYReader reader;
	// 把路径改为自己存放文件的路径
	reader.read ("monkey.ply", *cloud); // 记住要事先下载这个数据集！
	std::cerr << "PointCloud before filtering: " << cloud->width * cloud->height 
		 << " data points (" << pcl::getFieldsList (*cloud) << ").";
	// 创建滤波器对象
	pcl::VoxelGrid<pcl::PCLPointCloud2> sor;
	sor.setInputCloud (cloud);
	sor.setLeafSize (0.1f, 0.1f, 0.1f);
	sor.filter (*cloud_filtered);
	std::cerr << "PointCloud after filtering: " << cloud_filtered->width * cloud_filtered->height 
     << " data points (" << pcl::getFieldsList (*cloud_filtered) << ").";
	pcl::PCDWriter writer;
	writer.write("monkey2.pcd", *cloud_filtered);
#endif

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_in(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_out(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_icp(new pcl::PointCloud<pcl::PointXYZ>);
	// load file
#if FILTER
	if (pcl::io::loadPCDFile("monkey2.pcd", *cloud_in) == -1) {
		cout << "could not find the file" << endl;
		return -1;
	}
#else
	if (pcl::io::loadPLYFile("monkey.ply", *cloud_in) == -1) {
		cout << "could not find the file" << endl;
		return -1;
	}
#endif
	// 可视化原始点云
	pcl::visualization::PCLVisualizer viewer("demo");
	int v1 = 0;
	int v2 = 1;
	viewer.createViewPort(0.0, 0.0, 0.5, 1.0, v1);
	viewer.createViewPort(0.5, 0.0, 1.0, 1.0, v2);
	// 设置初始点云
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ>cloud_in_color(cloud_in,
		255, 255, 255);
	viewer.addPointCloud(cloud_in, cloud_in_color, "cloud_in_v1", v1);
	viewer.addPointCloud(cloud_in, cloud_in_color, "cloud_in_v2", v2);
	// 定义齐次变化矩阵
	Eigen::Matrix4d transformation_matrix = Eigen::Matrix4d::Identity();
	double theta = M_PI / 20; // 设置旋转弧度为9°
	transformation_matrix(0, 0) = cos(theta);
	transformation_matrix(0, 1) = -sin(theta);
	transformation_matrix(1, 0) = sin(theta);
	transformation_matrix(1, 1) = cos(theta);
	// 设置平移向量
	transformation_matrix(0, 3) = 0.0;
	transformation_matrix(1, 3) = 0.0;
	transformation_matrix(2, 3) = 0.0;
	pcl::transformPointCloud(*cloud_in, *cloud_out, transformation_matrix);
	// 设置初始转换后的点云
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ>cloud_out_color(cloud_out,
		20, 180, 20);
	viewer.addPointCloud(cloud_out, cloud_out_color, "cloud_out_v1", v1);

	pcl::console::TicToc time;
	mypcl::ICP<pcl::PointCloud<pcl::PointXYZ>::Ptr, pcl::PointCloud<pcl::PointXYZ>::Ptr> icp(cloud_in,cloud_out);
	//icp.SetEpsilon(1e-5);
	icp.SetIterations(50);
	icp.SetCorrThreshold(5);
	icp.SetMethodMS(mypcl::SVD);
	icp.SetMethodCP(mypcl::KD_TREE);
	time.tic(); // 开始计时
	icp.align(cloud_icp);
	int iteration_time = icp.getIterations();
	float SME = icp.getSquareMeanError();
	cout << "size of cloud points : " << cloud_in->points.size() << endl;
	cout << "total iteration times : " << iteration_time << endl;
	cout << "final square mean error : " << SME << endl;
	//cout << "score:" << icp.getFitnessScore() << endl;
	cout << "align time : " << time.toc() << endl;
	// 设置配准后的点云
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ>cloud_ICP_color(cloud_icp,
		180, 20, 20);
	viewer.addPointCloud(cloud_icp, cloud_ICP_color, "cloud_ICP_v2", v2);

	
	while (!viewer.wasStopped()) {
		viewer.spinOnce();
	}

	return 0;
}