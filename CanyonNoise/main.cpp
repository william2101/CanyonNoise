// Program is designed to downsample LiDAR pointcloud data
// Data is obtained from: https://portal.opentopography.org/datasets

#include <iostream>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

//Point structure used to generate the template class PointCloud
struct Point
{
	Point() : x(0), y(0), z(0) {}
	Point(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
	void print(){ cout << setprecision(10) << "(" << x << "," << y << "," << z <<")" << endl;}
    static Point mins (const Point& p1, const Point& p2); // return a point object with the smallest x,y,z componets between 2 points
    static Point maxes (const Point& p1, const Point& p2); // return a point object with the largest x,y,z componets between 2 points

    //Point members
	double x, y, z; 
};

Point Point::mins (const Point& p1, const Point& p2) {
    Point p;

    p.x = (p1.x < p2.x ? p1.x : p2.x);
    p.y = (p1.y < p2.y ? p1.y : p2.y);
    p.z = (p1.z < p2.z ? p1.z : p2.z);

    return p;
}

Point Point::maxes (const Point& p1, const Point& p2) {
    Point p;

    p.x = (p1.x > p2.x ? p1.x : p2.x);
    p.y = (p1.y > p2.y ? p1.y : p2.y);
    p.z = (p1.z > p2.z ? p1.z : p2.z);

    return p;
}

template <class PointT> class PointCloud {
	public:
	PointCloud() {
		data = NULL;
		size = 0;
	}

	PointCloud( const int size ) {
		data = NULL;
		allocate(size);
	}

    //PointCloud constructor from an input csv file
    PointCloud(string fileName) {

        ifstream f;
        f.open(fileName);

        if(!f){
                cout << "Error in reading input" << endl;
                std::exit(EXIT_FAILURE);
        }

        //find # of entries
        string unused;
        int numLines = 0;
        while ( std::getline(f, unused) )
            ++numLines;

        //return to start of file
        f.clear();
        f.seekg (0, ios::beg);

        //Allocate data matrix
        //Assumes 1 header line for now
        allocate(numLines - 1);

        string line;
        string val;
        double a[3] = {0, 0, 0};

        double element;

        getline(f, line); // skip the first line

        int i = 0;
        while (getline (f, line)) {         /* read each line */

            stringstream s(line);

            //Scan through each value of line
            int j = 0;
            while (getline (s, val, ',')){   /* for each value */
                if (!val.empty() && j <= 2){
                    a[j] = stod(val, 0);
                    j++; 
                }
                Point p_tmp(a[0], a[1], a[2]);
                data[i] = p_tmp;
            }
            i++;
        }
        f.close();
    
    }

	~PointCloud() {	deallocate(); }

	void insert(PointT p, int i){ data[i] = p; }

	int getSize() { return size; }

	void print() {
		for (int i = 0; i < size; i++){
			data[i].print();
		}
	}

    void write(string fileName) {
        ofstream f;

        f.open(fileName);
        f << "x,y,z" << endl;
        for (int i = 0; i < size; i++){
            f << setprecision(10) << data[i].x << "," << data[i].y << "," << data[i].z << endl;
        }
        f.close();
    }

    Point minVals() {
        Point p = data[0];

        for (int i = 1; i < size; i++){
            p = Point::mins(p, data[i]);
        }

        return p;
    }

    Point maxVals() {
        Point p = data[0];

        for (int i = 1; i < size; i++){
            p = Point::maxes(p, data[i]);
        }

        return p;
    }

    void downsample( PointCloud<PointT>& ds, int gridSize){

        int count = 0;
        int j = 0;
        Point p_mins = ds.minVals();
        Point p_maxes = ds.maxVals();
        double delta_x = (p_maxes.x - p_mins.x)/gridSize;
        double delta_y = (p_maxes.y - p_mins.y)/gridSize;
        double avg_x = 0, avg_y = 0, avg_z = 0;
        bool check = false;

        for (double x_ = p_mins.x; x_ < p_maxes.x; x_ += delta_x){
            for (double y_ = p_mins.y; y_ < p_maxes.y; y_ += delta_y){
                for (int i = 0; i < size; i++){
                    if (ds.data[i].x >= x_ && ds.data[i].x <= x_ + delta_x && ds.data[i].y >= y_ && ds.data[i].y <= y_ + delta_y){
                        avg_x += ds.data[i].x;
                        avg_y += ds.data[i].y;
                        avg_z += ds.data[i].z;
                        count += 1;
                        check = true;
                    } 
                }
                if (check){
                    Point p_tmp(avg_x/count, avg_y/count, avg_z/count);
                    insert(p_tmp, j);
                    avg_x = 0, avg_y = 0, avg_z = 0;
                    j += 1;
                    count = 0;
                    check = false;
                }
            }
        }
        size = j-1;
    }

	private:
	void allocate( const int& size ) {
		// Remember dimensions
		this->size = size;

		// Allocate
		data = new PointT[size];
	}

	void deallocate()
	{
		if ( NULL == data ) {// Nothing to do
			return;
		}

		// Free the memory
		delete[] data;

		// Reset
		size = 0;
		data = NULL;
	}

	int size;
	PointT* data;

};

int main()
{
    PointCloud<Point> PC("test_points.txt");
    PointCloud<Point> PC_DS(PC.getSize()); //Downsampled point cloud, allocate with initial size of input point cloud

    //Second parameter sets grid size
    PC_DS.downsample(PC, 40);

    PC_DS.write("test_out.csv");

    cout << "Original size (in # of points): " << PC.getSize() << endl;
    cout << "Downsampled size (in # of points): " << PC_DS.getSize() << endl;

}