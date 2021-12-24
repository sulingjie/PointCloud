#include <iostream>
#include <vector>
#include "MyGlobalICP.h"
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>

using namespace std;

int main(int argc, char** argv) {
	vector<pcl::PointCloud<pcl::PointXYZ>> clouds;
	vector<pcl::PointCloud<pcl::PointXYZ>>clouds_out;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud0(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud3(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::io::loadPCDFile("cloud0.pcd", *cloud0);
	pcl::io::loadPCDFile("cloud1.pcd", *cloud1);
	pcl::io::loadPCDFile("cloud2.pcd", *cloud2);
	pcl::io::loadPCDFile("cloud3.pcd", *cloud3);
	clouds.push_back(*cloud0);
	clouds.push_back(*cloud1);
	clouds.push_back(*cloud2);
	clouds.push_back(*cloud3);
	mypcl::GloballyICP globalICP(clouds,int(0));
	globalICP.setIterationNum(300);
	globalICP.setCorrThreshold(8);
	globalICP.setAlpha((float)0.001);
	globalICP.align(clouds_out);
	// 可视化
	// 可视化
	pcl::visualization::PCLVisualizer viewer("viewer");
	int v1 = 0;
	int v2 = 1;
	viewer.createViewPort(0.0, 0.0, 0.5, 1.0, v1);
	viewer.createViewPort(0.5, 0.0, 1.0, 1.0, v2);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud0Color(cloud0, 255, 255, 255);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud1Color(cloud1, 0, 0, 255);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud2Color(cloud2, 255, 0, 0);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud3Color(cloud3, 0, 255, 0);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud0outColor(clouds_out[0].makeShared(), 255, 255, 255);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud1outColor(clouds_out[1].makeShared(), 0, 0, 255);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud2outColor(clouds_out[2].makeShared(), 255, 0, 0);
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud3outColor(clouds_out[3].makeShared(), 0, 255, 0);
	viewer.addPointCloud(cloud0, cloud0Color, "cloud0",v1);
	viewer.addPointCloud(cloud1, cloud1Color, "cloud1",v1);
	viewer.addPointCloud(cloud2, cloud2Color, "cloud2",v1);
	viewer.addPointCloud(cloud3, cloud3Color, "cloud3",v1);
	viewer.addPointCloud(clouds_out[0].makeShared(), cloud0outColor, "cloud0out", v2);
	viewer.addPointCloud(clouds_out[1].makeShared(), cloud1outColor, "cloud1out", v2);
	viewer.addPointCloud(clouds_out[2].makeShared(), cloud2outColor, "cloud2out", v2);
	viewer.addPointCloud(clouds_out[3].makeShared(), cloud3outColor, "cloud3out", v2);
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cloud0");
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cloud1");
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cloud2");
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cloud3");
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cloud0out");
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cloud1out");
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cloud2out");
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cloud3out");
	while (!viewer.wasStopped()) {
		viewer.spinOnce();
	}
	return 0;
}