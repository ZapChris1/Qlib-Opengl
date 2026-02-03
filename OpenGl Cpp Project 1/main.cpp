#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
//opens windows
#include "dependencies/include/glad/glad.h"
//drawing
#include <cmath>
//for math(sin cos etc)
#include <iostream>

#include <vector>
#include <algorithm>
//for std::cout
#include <cstring>
//for memcpy

//this is pre processing and checks to see if M_Pi exists
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 6.28318530717958647692
#endif // !M_PI_2


#define SIZE 11
#define VERTEX_SIZE 7
#define FRAME_RATE 60
#define DEFAULT_COLOR {0.8f, 0.3f, 0.02f, 1.0f}
float FRAME_RATE_DIVISOR = (1.0 / FRAME_RATE);
//#define NEAR_PLANE 0.1
//#define FAR_PLANE 100
//#define FOV 90


const char* vertexShaderSource = "#version 330 core\n"
"layout (location = 0) in vec4 aPos;\n"
"layout (location = 1) in vec4 in_Color;\n"
"out vec4 FragPosition;\n"
"out vec4 FragColor;\n"
"uniform mat4 u_perspectiveMatrix;\n"
"uniform mat4 u_viewMatrix;\n"
"uniform mat4 u_positionMatrix;\n"
"void main()\n"
"{\n"
"	gl_Position = u_perspectiveMatrix * u_viewMatrix *  u_positionMatrix * vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
"	FragPosition = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
"   FragColor = in_Color;\n"
"}\0";



//transparency compares previous color to the new one (doesnt do anything unless opaque shapes are behind)
const char* fragmentShaderSource = "#version 330 core\n"
"in vec4 FragPosition;\n"
"in vec4 FragColor;\n"
"out vec4 out_Color;\n"
"void main()\n"
"{\n"
"	out_Color = vec4(FragColor.x, FragColor.y, FragColor.z, FragColor.w);\n"
"}\n\0";
//FragColor = vec4(0.8f, 0.3f, 0.02f, 1.0f);\n"

//these are not methods since they may be used in the future for the sphere class
void circleVerticesGenerator(unsigned int pointAmount, unsigned int lineAmnt, GLfloat* holder, unsigned int startIndex, float radius, float centerX, float centerY, GLuint* indexList, int indiceIndex) {
    for (unsigned int i = 0; i < pointAmount - 1; i++) {
        float pointNum = (float)i / (pointAmount - 1) * 2 * M_PI;
        float x = cos(pointNum) * radius + centerX;
        float y = sin(pointNum) * radius + centerY;
        holder[((i - 0) * VERTEX_SIZE) + startIndex] = x;
        holder[((i - 0) * VERTEX_SIZE + 1) + startIndex] = y;
        holder[((i - 0) * VERTEX_SIZE + 2) + startIndex] = 0.0f;
        //indexList[indiceIndex + i];
        //indexList[indiceIndex + i];
    }
    float pointNum = 1 * 2 * M_PI;
    float x = cos(pointNum) * radius + centerX;
    float y = sin(pointNum) * radius + centerY;
    holder[((pointAmount - 1) * VERTEX_SIZE) + startIndex] = x;
    holder[((pointAmount - 1) * VERTEX_SIZE + 1) + startIndex] = y;
    holder[((pointAmount - 1) * VERTEX_SIZE + 2) + startIndex] = 0.0f;
    //indexList[indiceIndex + (lineAmnt - 1)];
}

void ringVerticesGenerator(unsigned int pointAmount, unsigned int lineAmnt, GLfloat* holder, unsigned int startIndex, float radius, float centerX, float centerY, GLuint* indexList, int indiceIndex) {
    for (unsigned int i = 0; i < pointAmount - 1; i++) {
        float pointNum = (float)i / (pointAmount - 1) * 2 * M_PI;
        float x = cos(pointNum) * radius + centerX;
        float y = sin(pointNum) * radius + centerY;
        holder[((i - 0) * VERTEX_SIZE) + startIndex] = x;
        holder[((i - 0) * VERTEX_SIZE + 1) + startIndex] = y;
        holder[((i - 0) * VERTEX_SIZE + 2) + startIndex] = 0.0f;
        //indexList[indiceIndex + i];
        //indexList[indiceIndex + i];
    }
    float pointNum = 1 * 2 * M_PI;
    float x = cos(pointNum) * radius + centerX;
    float y = sin(pointNum) * radius + centerY;
    holder[((pointAmount - 1) * VERTEX_SIZE) + startIndex] = x;
    holder[((pointAmount - 1) * VERTEX_SIZE + 1) + startIndex] = y;
    holder[((pointAmount - 1) * VERTEX_SIZE + 2) + startIndex] = 0.0f;
    //indexList[indiceIndex + (lineAmnt - 1)];
}

void filledCircleVerticesGenerator(unsigned int pointAmount, GLfloat* holder, unsigned int startIndex, float radius, float centerX, float centerY) {
    holder[0 + startIndex] = centerX;
    holder[1 + startIndex] = centerY;
    holder[2 + startIndex] = 0.0f;
    for (unsigned int i = 3; i < pointAmount - 1; i++) {
        float pointNum = (float)i / (pointAmount - 1) * 2 * M_PI;
        float x = cos(pointNum) * radius + centerX - 0.5f;
        float y = sin(pointNum) * radius + centerY - 0.5f;
        holder[((i - 0) * VERTEX_SIZE) + startIndex] = x;
        holder[((i - 0) * VERTEX_SIZE + 1) + startIndex] = y;
        holder[((i - 0) * VERTEX_SIZE + 2) + startIndex] = 0.0f;
    }
    float pointNum = 1 * 2 * M_PI;
    float x = cos(pointNum) * radius + centerX - 0.5f;
    float y = sin(pointNum) * radius + centerY - 0.5f;
    holder[((pointAmount - 1) * VERTEX_SIZE) + startIndex] = x;
    holder[((pointAmount - 1) * VERTEX_SIZE + 1) + startIndex] = y;
    holder[((pointAmount - 1) * VERTEX_SIZE + 2) + startIndex] = 0.0f;
}

class Vector3 {
public:
    float x, y, z;
    void flip() {
        x = -x;
        y = -y;
        z = -z;
    }
    Vector3 operator*(const float w) {
        Vector3 outputVector;
        outputVector.x = x * w;
        outputVector.y = y * w;
        outputVector.z = z * w;
        return outputVector;
    }
};

struct Vector4 {
    float x, y, z, w;
};

class matrix3x3 {
public:
    double m[3][3];
    Vector3 operator*(const Vector3& v) {
        Vector3 newVector;
        newVector.x = v.x * m[0][0] + v.y * m[0][1] + v.z * m[0][2];
        newVector.y = v.x * m[1][0] + v.y * m[1][1] + v.z * m[1][2];
        newVector.z = v.x * m[2][0] + v.y * m[2][1] + v.z * m[2][2];
        return newVector;
    }
};
//naive unsigned integer exponents
unsigned int pow(unsigned int base, unsigned int exponent) {
    unsigned int result = 1;
    for (unsigned int i = 0; i < exponent; i++) {
        result *= base;
    }
    return result;
}
//naive signed integer exponents
int pow(int base, int exponent) {
    int result = 1;
    for (int i = 0; i < exponent; i++) {
        result *= base;
    }
    return result;
}

class Quaternion {
public:
    Quaternion(double W, double X, double Y, double Z) {
        w = W;
        x = X;
        y = Y;
        z = Z;
    }

    Quaternion() {
        w = 1.0;
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }

    //Hamilton Product
    Quaternion operator*(const Quaternion& q2) {
        return {
            w * q2.w - x * q2.x - y * q2.y - z * q2.z, //w
                w * q2.x + x * q2.w + y * q2.z - z * q2.y, //x
                w * q2.y - x * q2.z + y * q2.w + z * q2.x, //y
                w * q2.z + x * q2.y - y * q2.x + z * q2.w //z
        };
    }
    Vector3 getForward() {
        Vector3 outputVector;
        outputVector.x = -2 * x * z + 2 * y * w;
        outputVector.y = -2 * y * z - 2 * x * w;
        outputVector.z = -1 + 2 * x * x + 2 * y * y;
        return outputVector;
    }
    Vector3 getRight() {
        Vector3 outputVector;
        outputVector.x = 1 - 2 * y * y - 2 * z * z;
        outputVector.y = 2 * x * y - 2 * z * w;
        outputVector.z = 2 * x * z + 2 * y * w;
        return outputVector;
    }
    Vector3 getUp() {
        Vector3 outputVector;
        outputVector.x = 2 * x * y + 2 * z * w;
        outputVector.y = 1 - 2 * x * x - 2 * z * z;
        outputVector.z = 2 * y * z - 2 * x * w;
        return outputVector;
    }
    //Scalar
    Quaternion operator*(const double k) {
        return {
            (w * k),
            (x * k),
            (y * k),
            (z * k)
        };
    }
    //real dotproduct
    Quaternion operator+(const Quaternion& q2) {
        return {
            w + q2.w,
            x + q2.x,
            y + q2.y,
            z + q2.z
        };
    }
    Quaternion operator-(const Quaternion& q2) {
        return {
            w - q2.w,
            x - q2.x,
            y - q2.y,
            z - q2.z
        };
    }
    //the inverse for 3D purposes
    Quaternion conjugate() {
        return {
            w,
            x *= -1,
            y *= -1,
            z *= -1
        };
    }

    static double dotProduct(const Quaternion& q1, const Quaternion& q2) {
        //scalar returned: + not , on newlines
        return {
           q1.w * q2.w +
           q1.x * q2.x +
           q1.y * q2.y +
           q1.z * q2.z
        };
    }

    static double magnitiude(const Quaternion& q1) {
        return {
            std::sqrt((q1.w * q1.w) +
            (q1.x * q1.x) +
            (q1.y * q1.y) +
            (q1.z * q1.z))
        };
    }

    void normalize() {
        double magnitude =
            std::sqrt((w * w) +
                (x * x) +
                (y * y) +
                (z * z));
        w = (w / magnitude);
        x = (x / magnitude);
        y = (y / magnitude);
        z = (z / magnitude);
    }

    matrix3x3 getMatrix() {
        matrix3x3 m;
        m.m[0][0] = 1 - 2 * y * y - 2 * z * z;
        m.m[0][1] = 2 * x * y - 2 * z * w;
        m.m[0][2] = 2 * x * z + 2 * y * w;

        m.m[1][0] = 2 * x * y + 2 * z * w;
        m.m[1][1] = 1 - 2 * x * x - 2 * z * z;
        m.m[1][2] = 2 * y * z - 2 * x * w;


        m.m[2][0] = 2 * x * z - 2 * y * w;
        m.m[2][1] = 2 * y * z + 2 * x * w;
        m.m[2][2] = 1 - 2 * x * x - 2 * y * y;
        return m;
        //struct is a wrapper to return the double[3][3] as a nonpointer so it stays on the stack (faster)
    }

    double w, x, y, z;
};

class SLERP {
public:
    SLERP(Quaternion quaternion1, Quaternion quaternion2) {
        q1 = quaternion1;
        q2 = quaternion2;
        time = 0.0;
        calcTheta();
    }

    void reboot(Quaternion quaternion1, Quaternion quaternion2) {
        q1 = quaternion1;
        q2 = quaternion2;
        time = 0.0;
        calcTheta();
    }

    Quaternion getQ3(double t) {
        return calcQ3(t);
    }

    void setTime(double t) {
        time = t;
    }
    
    void addTime(double t) {
        time += t;
        if (t > 1) {
            t = t - 1;
        }
    }

    Quaternion getQ3(void) {
        return calcQ3(time);
    }

    inline double getTheta(void) {
        return theta;
    }
    inline Quaternion getQ1(void) {
        return q1;
    }
    inline Quaternion getQ2(void) {
        return q2;
    }
private:
    void calcTheta(void) {
        double dotProduct = Quaternion::dotProduct(q1, q2);
        if (dotProduct > 1.0) dotProduct = 1.0;
        if (dotProduct < -1.0) dotProduct = -1.0;
        if (dotProduct == 0.0) {
            theta = acos(Quaternion::dotProduct(q1, q2 * -1));
            sinTheta = sin(theta);
        }
        else {
            theta = acos(dotProduct);
            sinTheta = sin(theta);
        }
    }

    Quaternion calcQ3(double t) {
        if (theta == 0.0) {
            Quaternion q3 = (q1 * (1 - t)) + (q2 * t);
            return q3;
        }
        double scalar1 = sin((1 - t) * theta) / sinTheta;
        double scalar2 = sin(t * theta) / sinTheta;
        Quaternion scaledQ1 = q1 * scalar1;
        Quaternion scaledQ2 = q2 * scalar2;
        Quaternion q3 = scaledQ1 + scaledQ2;
        return q3;
    }
    double theta;
    Quaternion q1;
    Quaternion q2;
    double sinTheta;
    double time;
};

//put 0 for indexLength to just default to what was previously there
class Shape {
public:
    //length is amt of vertices and indexLength is amt of indices
    Shape(GLenum input_draw_type, unsigned int length, unsigned int indexLength, float X = 0, float Y = 0, Vector4 Color = DEFAULT_COLOR) {
        ready = false;
        drawType = input_draw_type;
        x = X;
        y = Y;
        color = Color;
        ++shapeCount;
        if (!indexLength) {
            //fires
            this->len = length;
        }
        else {
            this->len = indexLength;
        }
        length *= VERTEX_SIZE;
        this->length = length;
        //std::cout << "parentListSize: " << parentListSize << "  length: " << length << std::endl;
        //std::cout << "indexListSize: " << indexListSize << "  len: " << len << std::endl;
        //std::cout << "shapeCount: " << shapeCount << std::endl;
        adjustorValue = length;
        lastAdjustorValue = adjustorValue;
        currentIndex = parentListSize;
        shapeIndex = shapeCount - 1;
        indiceIndex = indexListSize;
        if (parentList != nullptr) {
            GLfloat* newList = new GLfloat[parentListSize + length];
            memcpy(newList, parentList, parentListSize * sizeof(GLfloat));

            delete[] parentList;

            parentList = newList;
            newList = nullptr;

            parentListSize += length;
        }
        else {
            parentList = new GLfloat[parentListSize + length];
            parentListSize += length;
        }
        if (indexList != nullptr) {
            GLuint* newIndexList = new GLuint[indexListSize + len];
            memcpy(newIndexList, indexList, indexListSize * sizeof(GLuint));

            delete[] indexList;

            indexList = newIndexList;
            newIndexList = nullptr;

            indexListSize += len;
        }
        else {
            indexList = new GLuint[indexListSize + len];
            indexListSize += len;
        }

        if (shapeList != nullptr) {
            Shape** newList = new Shape * [shapeCount];
            memcpy(newList, shapeList, (shapeCount - 1) * sizeof(Shape*));  // COPY ONLY shapeCount-1
            delete[] shapeList;
            shapeList = newList;
            newList = nullptr;
            //shapecount was handled prior
        }
        else {
            shapeList = new Shape * [shapeCount];
        }
        shapeList[shapeCount - 1] = this;  // Insert this instance
        //std::cout << "End Create" << std::endl;
        std::cout << "NewSize: " << parentListSize << std::endl;
    }

    static void setReady(bool boolean) {
        ready = boolean;
    }

    static GLfloat* getParentList(void) {
        return parentList;
    }
    static GLuint* getIndexList(void) {
        return indexList;
    }

    static unsigned int getParentListSize(void) {
        return parentListSize;
    }

    static unsigned int getIndexListSize(void) {
        return indexListSize;
    }

    void setVertice(unsigned int verticeIndex, GLfloat verticeData) {
        if (verticeIndex < length) parentList[verticeIndex + currentIndex] = verticeData;
        else std::cout << "something odd happened line: " << __LINE__ << " of file: " << __FILE__ << std::endl;
    }

    void drawSelf(void) {
        if (drawType) {
            if (ready) {
                glDrawArrays(drawType, currentIndex / VERTEX_SIZE, length / VERTEX_SIZE);
            }
        }
        else std::cout << "something odd happened line: " << __LINE__ << " of file: " << __FILE__ << std::endl;
    }
    //chaningpoint is the currentindex of the element removed

    float X(void) {
        return x;
    }
    float Y(void) {
        return y;
    }


    void destroy(void) {
        this->~Shape();
        destroyed = true;
    }

    ~Shape(void) {
        if (!destroyed) {
            std::cout << "Start Destroy" << std::endl;
            --shapeCount;
            adjustorValue -= length;
            if (shapeCount < 1) {
                std::cout << "Destroy Last" << std::endl;
                std::cout << "parentListSize: " << parentListSize << "  length: " << length << std::endl;
                std::cout << "indexListSize: " << indexListSize << "  len: " << len << std::endl;
                std::cout << "shapeCount: " << shapeCount << std::endl;
                //deallocate memory
                delete[] parentList;
                parentList = nullptr;
                delete[] indexList;
                indexList = nullptr;
                delete[] shapeList;
                shapeList = nullptr;
                std::cout << "Destroy Completed" << std::endl;
            }
            else if (parentList != nullptr || shapeList != nullptr) {
                if (parentList != nullptr) {
                    std::cout << "Destroy Left" << std::endl;
                    std::cout << "parentListSize: " << parentListSize << "  length: " << length << std::endl;
                    std::cout << "indexListSize: " << indexListSize << "  len: " << len << std::endl;
                    std::cout << "shapeCount: " << shapeCount << std::endl;
                    //
                    GLfloat* newList = new GLfloat[parentListSize - length];
                    GLuint* newIndexList = new GLuint[indexListSize - len];
                    std::cout << "Destroy Left 1.1" << std::endl;
                    //copy first part before the removed part
                    //if currentIndex != 0
                    if (currentIndex) memcpy(newList, parentList, sizeof(GLfloat) * currentIndex);
                    std::cout << "Destroy Left 1.2" << std::endl;
                    //copy second part past removed part
                    memcpy(newList, parentList + ((currentIndex + length) - 1), sizeof(GLfloat) * (parentListSize - length));
                    std::cout << "Destroy Left 1.3" << std::endl;
                    //copy first part then manually copy rest (done manually so that I can subtract by len)
                    //if indiceIndex != 0
                    if (indiceIndex) memcpy(newIndexList, indexList, sizeof(GLuint) * indiceIndex);
                    std::cout << "Destroy Left 1.6" << std::endl;
                    //temporary value made for optimization

                    int tempValue = ((indiceIndex + len) - 0);
                    std::cout << "Destroy Left 1.8" << std::endl;
                    //
                    for (int i = indiceIndex + len; i < indexListSize; ++i) {
                        //std::cout << (indexListSize - (i + 1)) << std::endl;
                        newIndexList[i] = indexList[i] - len;
                    }
                    std::cout << "Destroy Left end" << std::endl;
                    //cleanup
                    delete[] parentList;
                    parentList = newList;
                    delete[] indexList;
                    indexList = newIndexList;
                    parentListSize -= length;
                    indexListSize -= len;

                }
                if (shapeList != nullptr) {
                    //we subtract shapecount by 1 less since its alr -- earlier
                    Shape** newList = new Shape * [shapeCount + 1];  // already --shapeCount at top
                    for (int i = 0, j = 0; i < shapeCount + 1; i++) {
                        if (shapeList[i] == this) continue;
                        if (i > shapeIndex) {
                            //adjust the values of each shape
                            shapeList[i]->currentIndex -= length;
                            shapeList[i]->shapeIndex -= 1;
                            shapeList[i]->indiceIndex -= len;
                        }
                        newList[j++] = shapeList[i];
                    }
                    delete[] shapeList;
                    shapeList = newList;
                }
            }
            ready = false;
            std::cout << "End Destroy" << std::endl;
        }
    }
    //RGB then opacity (1-0 for all fields
    Vector4 color;
protected:
    bool destroyed = false;
    static bool ready;
    static unsigned int shapeCount;
    static signed int adjustorValue;
    static GLfloat* parentList;
    static unsigned int parentListSize;
    static GLuint* indexList;
    static GLfloat* adjustedindexList;
    static unsigned int indexListSize;
    static Shape** shapeList;

    float x;
    float y;
    float rotation = 0;
    unsigned int currentIndex;
    unsigned int shapeIndex;
    unsigned int indiceIndex;
    unsigned int length;
    unsigned int len;
    signed int lastAdjustorValue;
    GLenum drawType;
};
bool Shape::ready = true;
unsigned int Shape::shapeCount = 0;
signed int Shape::adjustorValue = 0;
GLfloat* Shape::parentList = nullptr;
GLuint* Shape::indexList = nullptr;
GLfloat* Shape::adjustedindexList = nullptr;
unsigned int Shape::parentListSize = 0;
unsigned int Shape::indexListSize = 0;
Shape** Shape::shapeList = nullptr;

//base for 2d circles
class Circle : public Shape {
public:
    Circle(GLenum input_draw_type, unsigned int len, float input_radius, float X, float Y, bool input_fill) : Shape(input_draw_type, len, 0, X, Y) {
        std::cout << len << std::endl;
        filled = input_fill;
        radius = input_radius;
        //drawCircle();
    }

    void drawCircle(void) {
        circleVerticesGenerator(length, len, parentList, currentIndex, radius, x, y, indexList, indiceIndex);
        //if (filled) filledCircleVerticesGenerator(length / VERTEX_SIZE, parentList, currentIndex, radius, x, y);
        //else circleVerticesGenerator(length / VERTEX_SIZE, parentList, currentIndex, radius, x, y);   
    }

    void setRadius(float new_radius) {
        radius = new_radius;
        drawCircle();
    }

    void setXY(float X, float Y) {
        x = X;
        y = Y;
        drawCircle();
    }

    void setX(float X) {
        x = X;
        drawCircle();
    }

    void setY(float Y) {
        y = Y;
        drawCircle();
    }

protected:
    float radius;
    bool filled;
};

//base for 2d shapes with rotation
class Polygon : public Shape {
public:
    Polygon(GLenum input_draw_type, float len, float X, float Y) : Shape(input_draw_type, len, 0, X, Y) {
        localList = new GLfloat[length];
    }

    void copyToLocalList(void) {
        //must be done only once and after drawing the shape
        memcpy(localList, parentList + currentIndex, length * sizeof(GLfloat));
        localListInitialized = true;
    }

    void setRotation(float Rotation) {
        if (localListInitialized) {
            rotation = Rotation;
            rotatePoints2D();
        }
    }

    void copyToParentList(void) {
        memcpy(parentList + currentIndex, localList, length * sizeof(GLfloat));
    }

    void changeRotation(float deltaRotation) {
        if (localListInitialized) {
            rotation += deltaRotation;
            rotatePoints2D();
        }
    }
    float returnRotation(void) {
        return rotation;
    }

    void rotatePoints2D() {
        for (unsigned int index = 0; index < length / VERTEX_SIZE; index++) {
            float& dx = localList[(index * VERTEX_SIZE)];
            float& dy = localList[(index * VERTEX_SIZE) + 1];
            const float Cos = cos(rotation);
            const float Sin = sin(rotation);
            float fx = dx * Cos - dy * Sin;
            float fy = dx * Sin + dy * Cos;
            parentList[(index * VERTEX_SIZE) + currentIndex] = fx + x;
            parentList[(index * VERTEX_SIZE) + 1 + currentIndex] = fy + y;
        }
    }


    ~Polygon(void) {
        delete[] localList;
    }
protected:
    float rotation = 0;
    bool localListInitialized = false;
    GLfloat* localList;
};

//base for 3d shapes with rotation
class Polygon3D : public Shape {
public:
    Polygon3D(GLenum input_draw_type, int length, int indexLength, float X, float Y, float Z, Vector4 Color) : Shape(input_draw_type, length, indexLength, X, Y, Color) {
        localList = new GLfloat[length * VERTEX_SIZE];
        z = Z;
    }

    //initializes the localList
    void copyToLocalList(void) {
        ;
        memcpy(localList, parentList + currentIndex, length * sizeof(GLfloat));
        //localListInitialized = true;
    }

    //simply uses the quaternion stored
    void updateRotation() {
        if (true) {
            //take from the localList that preserves offset
            //rotatePoints3D(parentList, localList, currentIndex, length / VERTEX_SIZE, rotationQuaternion, x, y, z);
            rotatePoints3D(parentList, localList, currentIndex, length / VERTEX_SIZE, rotationQuaternion, x, y, z);
        }
    }

    void copyToParentList(void) {
        memcpy(parentList + currentIndex, localList, length * sizeof(GLfloat));
    }
    float Z(void) {
        return z;
    }
    ~Polygon3D(void) {
        delete[] localList;
    }

    //quaternion is made public to rotate and set via direct control for versatility
    Quaternion rotationQuaternion;
protected:

    void rotatePoints3D(GLfloat* holder, GLfloat* originalList, unsigned int startIndex, unsigned int pointAmount, Quaternion q1, float centerX, float centerY, float centerZ) {
        matrix3x3 m = q1.getMatrix();
        for (unsigned int index = 0; index < length / VERTEX_SIZE; index++) {
            float tempX = localList[(index * VERTEX_SIZE)] - x;
            float tempY = localList[(index * VERTEX_SIZE) + 1] - y;
            float tempZ = localList[(index * VERTEX_SIZE) + 2] - z;

            parentList[(index * VERTEX_SIZE) + currentIndex] = ((m.m[0][0] * tempX) + (m.m[0][1] * tempY) + (m.m[0][2] * tempZ)) + x;
            parentList[(index * VERTEX_SIZE) + 1 + currentIndex] = ((m.m[1][0] * tempX) + (m.m[1][1] * tempY) + (m.m[1][2] * tempZ)) + y;
            parentList[(index * VERTEX_SIZE) + 2 + currentIndex] = ((m.m[2][0] * tempX) + (m.m[2][1] * tempY) + (m.m[2][2] * tempZ)) + z;
        }
    }

    float rotation = 0;
    float z;
    bool localListInitialized = false;
    GLfloat* localList;
};

//base for 2d rectangles
class Rectangle : public Polygon {
public:
    Rectangle(GLenum input_draw_type, float X, float Y, float Height, float Width, float Depth) : Polygon(input_draw_type, 4, X, Y) {
        height = Height;
        width = Width;
        depth = Depth;
        drawRectangle();
    }

    void drawRectangle() {
        //this check is included here because unlike with quaternions its much cheaper to check and this is common enough to be considered a worthwhile usecase
        if (rotation == 0) {
            rectangleVerticesGenerator(parentList, currentIndex);
            copyToLocalList();
        }
        else {
            rectangleVerticesGenerator(parentList, currentIndex);
            rectangleVerticesGenerator(localList, 0);
            localListInitialized = true;
        }
    }

    void setXY(float X, float Y) {
        x = X;
        y = Y;
        drawRectangle();
    }

    void setX(float X) {
        x = X;
        drawRectangle();
    }

    void setY(float Y) {
        y = Y;
        drawRectangle();
    }

protected:
    //holder is the destination
    void rectangleVerticesGenerator(GLfloat* holder, unsigned int startIndex) {
        for (unsigned int i = 0; i < 4; i++) {
            //X and Y are the temporary values for the calculation while x and y are the properties of the class
            float X = (i < 2) ? x + (width / 2) : x - (width / 2);
            float Y = (i < 1 || i > 2) ? y + (height / 2) : y - (height / 2);
            holder[(i * VERTEX_SIZE) + 0 + startIndex] = X;
            holder[(i * VERTEX_SIZE) + 1 + startIndex] = Y;
            holder[(i * VERTEX_SIZE) + 2 + startIndex] = 0.0f;
            indexList[indiceIndex + i];
        }
    }

    float height;
    float width;
    float depth;
};

//base for 3d rectangles
class Rectangle3D : public Polygon3D {
public:
    Rectangle3D(GLenum input_draw_type, float X, float Y, float Z, float Height, float Width, float Depth, Vector4 Color = DEFAULT_COLOR) : Polygon3D(input_draw_type, 8, 24, X, Y, Z, Color) {
        height = Height;
        width = Width;
        depth = Depth;

        //called manually to prevent need for initialization check during runtime
        rectangleVerticesGenerator3D(parentList, currentIndex);
        //localList is a backup to preserve unrotated position (this preserves potential loss of data from floating point errors)
        copyToLocalList();
        //since order doesnt change we only need to generate indices once (if its needed in the future it will be made public for that usecase)
        rectangle3dIndexGenerator();
        localListInitialized = true;
        //drawRectangle();
    }

    //this function is assumed to run after localist has been initialized by copying from the parent list
    void drawRectangle() {
        //generate seperate list of vertices for each
        rectangleVerticesGenerator3D(parentList, currentIndex);
        //localList is a backup to preserve unrotated position (this preserves potential loss of data from floating point errors)
        rectangleVerticesGenerator3D(localList, 0);
        updateRotation();

    }

    void setXYZ(float X, float Y, float Z) {
        x = X;
        y = Y;
        z = Z;
        drawRectangle();
    }
    void setX(float X) {
        x = X;
        drawRectangle();
    }
    void setY(float Y) {
        y = Y;
        drawRectangle();
    }
    void setZ(float Z) {
        z = Z;
        drawRectangle();
    }

    //color is public to be modified the same way that the rotation is
    void update(void) {
        drawRectangle();
    }

    void setPositionAndColor(float X, float Y, float Z, Vector4 Color) {
        x = X;
        y = Y;
        z = Z;
        color = Color;
        drawRectangle();
    }

    //this version uses the indexBuffer 
    void drawSelfIndex(void) {
        if (drawType) {
            if (ready) {
                glDrawElements(GL_LINES, len, GL_UNSIGNED_INT, indexList + indiceIndex);
            }
        }
        else std::cout << "something odd happened line: " << __LINE__ << " of file: " << __FILE__ << std::endl;
    }
    void drawSelfIndex(bool iboBounded) {
        if (drawType) {
            if (!iboBounded) {
                glDrawElements(GL_LINES, len, GL_UNSIGNED_INT, indexList + indiceIndex);
            }
            else {
                //glDrawElements(GL_LINES, len, GL_UNSIGNED_INT, (void*)(indiceIndex * sizeof(GLuint)));
                glDrawElements(GL_LINES, len, GL_UNSIGNED_INT, (void*)(indiceIndex * sizeof(GLuint)));
            }
        }
        else std::cout << "something odd happened line: " << __LINE__ << " of file: " << __FILE__ << std::endl;
    }
    //chaningpoint is the currentindex of the element removed

protected:
    //holder is the destination
    //everything was hardcoded into it as constants for runtime efficency (its faster to read than generate new values)
    void rectangleVerticesGenerator3D(GLfloat* holder, unsigned int startIndex) {
        //x values
        //holder[0] = 1;
        holder[(0 * VERTEX_SIZE) + 0 + startIndex] = x + (width / 2);
        holder[(1 * VERTEX_SIZE) + 0 + startIndex] = x + (width / 2);
        holder[(2 * VERTEX_SIZE) + 0 + startIndex] = x - (width / 2);
        holder[(3 * VERTEX_SIZE) + 0 + startIndex] = x - (width / 2);
        holder[(4 * VERTEX_SIZE) + 0 + startIndex] = x + (width / 2);
        holder[(5 * VERTEX_SIZE) + 0 + startIndex] = x + (width / 2);
        holder[(6 * VERTEX_SIZE) + 0 + startIndex] = x - (width / 2);
        holder[(7 * VERTEX_SIZE) + 0 + startIndex] = x - (width / 2);
        //y values
        holder[(0 * VERTEX_SIZE) + 1 + startIndex] = y - (height / 2);
        holder[(1 * VERTEX_SIZE) + 1 + startIndex] = y + (height / 2);
        holder[(2 * VERTEX_SIZE) + 1 + startIndex] = y - (height / 2);
        holder[(3 * VERTEX_SIZE) + 1 + startIndex] = y + (height / 2);
        holder[(4 * VERTEX_SIZE) + 1 + startIndex] = y - (height / 2);
        holder[(5 * VERTEX_SIZE) + 1 + startIndex] = y + (height / 2);
        holder[(6 * VERTEX_SIZE) + 1 + startIndex] = y - (height / 2);
        holder[(7 * VERTEX_SIZE) + 1 + startIndex] = y + (height / 2);
        //z values
        holder[(0 * VERTEX_SIZE) + 2 + startIndex] = z - (depth / 2);
        holder[(1 * VERTEX_SIZE) + 2 + startIndex] = z - (depth / 2);
        holder[(2 * VERTEX_SIZE) + 2 + startIndex] = z - (depth / 2);
        holder[(3 * VERTEX_SIZE) + 2 + startIndex] = z - (depth / 2);
        holder[(4 * VERTEX_SIZE) + 2 + startIndex] = z + (depth / 2);
        holder[(5 * VERTEX_SIZE) + 2 + startIndex] = z + (depth / 2);
        holder[(6 * VERTEX_SIZE) + 2 + startIndex] = z + (depth / 2);
        holder[(7 * VERTEX_SIZE) + 2 + startIndex] = z + (depth / 2);
        //r values
        holder[(0 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(1 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(2 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(3 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(4 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(5 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(6 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(7 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        //g values
        holder[(0 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(1 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(2 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(3 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(4 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(5 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(6 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(7 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        //b values
        holder[(0 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(1 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(2 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(3 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(4 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(5 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(6 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(7 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        //opacity values
        holder[(0 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(1 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(2 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(3 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(4 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(5 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(6 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(7 * VERTEX_SIZE) + 6 + startIndex] = color.w;

    }
    void rectangle3dIndexGenerator(void) {
        //series of lines tells the GPU the order to connect points so we dont have to send duplicate vertices
        unsigned int currentIndexAdjusted = currentIndex / VERTEX_SIZE;
        GLuint tempIndexList[24] = {
            5 + currentIndexAdjusted, 7 + currentIndexAdjusted,
            7 + currentIndexAdjusted, 3 + currentIndexAdjusted,
            3 + currentIndexAdjusted, 1 + currentIndexAdjusted,
            3 + currentIndexAdjusted, 2 + currentIndexAdjusted,
            2 + currentIndexAdjusted, 0 + currentIndexAdjusted,
            2 + currentIndexAdjusted, 6 + currentIndexAdjusted,
            6 + currentIndexAdjusted, 7 + currentIndexAdjusted,
            6 + currentIndexAdjusted, 4 + currentIndexAdjusted,
            4 + currentIndexAdjusted, 5 + currentIndexAdjusted,
            4 + currentIndexAdjusted, 0 + currentIndexAdjusted,
            0 + currentIndexAdjusted, 1 + currentIndexAdjusted,
            1 + currentIndexAdjusted, 5 + currentIndexAdjusted };
        memcpy(indexList + indiceIndex, tempIndexList, sizeof(GLuint) * 24);
    }
    float height;
    float width;
    float depth;
    bool filled;
};



//make a sphere through the SphereConstructor class
class Sphere : public Polygon3D {
public:
    void setRadius(float new_radius) {
        radius = new_radius;
        drawSphere();
    }

    void setXYZ(float X, float Y, float Z) {
        offsetArray((X - x), (Y - y), (Z - z), parentList, currentIndex, length);
        offsetArray((X - x), (Y - y), (Z - z), localList, 0, length);
        x = X;
        y = Y;
        z = Z;
        //drawSphere();
    }

    void setX(float X) {
        offsetArray((X - x), 0, 0, parentList, currentIndex, length);
        offsetArray((X - x), 0, 0, localList, 0, length);
        x = X;
        //drawSphere();
    }

    void setY(float Y) {
        offsetArray(0, (Y - y), 0, parentList, currentIndex, length);
        offsetArray(0, (Y - y), 0, localList, 0, length);
        y = Y;
        //drawSphere();
    }

    void setZ(float Z) {
        offsetArray(0, 0, (Z - z), parentList, currentIndex, length);
        offsetArray(0, 0, (Z - z), localList, 0, length);
        z = Z;
        //drawSphere();
    }
    
    void drawSphere() {
        
        sphereVerticesGeneratorV2(parentList, currentIndex);
        //localList is the vertices unadjusted for position
        sphereVerticesGeneratorV2(localList, 0);
        updateRotation();
    }

    void drawSelfIndex(bool iboBounded) {
        if (drawType) {
            if (!iboBounded) {
                glDrawElements(GL_TRIANGLE_STRIP, len, GL_UNSIGNED_INT, indexList + indiceIndex);
            }
            else {
                //glDrawElements(GL_LINES, len, GL_UNSIGNED_INT, (void*)(indiceIndex * sizeof(GLuint)));
                glDrawElements(GL_TRIANGLE_STRIP, len, GL_UNSIGNED_INT, (void*)(indiceIndex * sizeof(GLuint)));
            }
        }
        else std::cout << "something odd happened line: " << __LINE__ << " of file: " << __FILE__ << std::endl;
    }
    
    ~Sphere() {
        delete[] verticeRingAmounts;
        delete[] uniqueValues;
    }

protected:
    friend class SphereConstructor;

    void offsetArray(float X, float Y, float Z, float* array, unsigned int currentIndex, unsigned int verticeAmount) {
        for (unsigned int i = 0; i < verticeAmount; i++) {
            unsigned int newIndex = (i * 3) + currentIndex;
            array[newIndex] += X;
            array[newIndex + 1] += Y;
            array[newIndex + 2] += Z;
        }
    }
    //like a UV sphere but with 4 more vertices on each ring starting with 1->4->8->12 etc it goes up to the middle ring then goes back down
    Sphere(GLenum input_draw_type, unsigned int length, unsigned int indexLength, unsigned int ringAmount, unsigned int* ringAmounts, unsigned int halfRingAmount, float input_radius, float X, float Y, float Z, Vector4 Color) : Polygon3D(input_draw_type, length, indexLength/*index*/, X, Y, Z, Color) {
        radius = input_radius;
        verticeRingAmounts = ringAmounts;
        this->ringAmount = ringAmount;
        this->halfRingAmount = halfRingAmount;
        uniqueValues = new float[ringAmount];
        uniqueValueGenerator(uniqueValues, ringAmount, halfRingAmount);
        sphereVerticesGeneratorV2(parentList, currentIndex);
        sphereIndicesGeneratorV3(); //safe but inefficent (overallocating memory)
        incrementer = 0;
    }

    //t is the smaller ring, z is the bigger ring, and f is the next ring, set f to true if there is a next ring and false if not //currentIndex is for the SPECIFIC ring not the currentIndex for the whole
    unsigned int fillRing(unsigned int t, unsigned int z, unsigned int t_total, unsigned int z_total, bool f, unsigned int currentIndex, bool next_ring_is_larger) { //currentIndex is for the SPECIFIC ring not the currentIndex for the whole
        indexList[currentIndex] = z; indexList[currentIndex + 1] = t; //initial 2 points for base of the first triangle
        currentIndex += 2;
        unsigned int t_counter = 4; unsigned int z_counter = 4;
        unsigned int corner_point_amount = (2 * (((z_total / 4) - 1) - 2)) + 3; //amount of triangles between degens
        indexList[currentIndex] = z + z_counter; currentIndex++; //initial 2 points for 1st ring
        indexList[currentIndex] = t + t_counter; currentIndex++;
        bool z_or_t = 0; //0 = z 1 = t
        unsigned int degenCounter = 0; //start on zero since technically first triangle hasnt been drawn yet
        unsigned int degenCounterMax = 9999; //logic for setting this comes later
        for (unsigned int i = 0; i < corner_point_amount; i++) { //first ring [1:]
            if (false) {
                indexList[currentIndex] = indexList[currentIndex - 2];
                t_counter = indexList[currentIndex - 2] - t;
            }
            else {
                if (z_or_t) { //t
                    t_counter += 8;
                    indexList[currentIndex] = t + t_counter;
                }
                else { //z
                    z_counter += 8;
                    indexList[currentIndex] = z + z_counter;
                }
            }
            currentIndex++;
            z_or_t = !z_or_t;
        }


        indexList[currentIndex] = indexList[currentIndex - 2];  //First Degenerate triangle
        t_counter = indexList[currentIndex - 2] - t;  z_or_t = !z_or_t;
        currentIndex++;
        z_counter += -6; t_counter += -6;
        if (z_or_t) {
            indexList[currentIndex] = t + t_counter; currentIndex++;
            indexList[currentIndex] = z + z_counter; currentIndex++; //initial 2 points for 2nd ring (4)
        }
        else { 
            indexList[currentIndex] = z + z_counter; currentIndex++; //initial 2 points for 2nd ring (4)
            indexList[currentIndex] = t + t_counter; currentIndex++;
        }
        for (unsigned int i = 0; i < corner_point_amount - 2; i++) { //2nd ring [1:-1] 
            if (false) {
                indexList[currentIndex] = indexList[currentIndex - 2];
                t_counter = indexList[currentIndex - 2] - t;
            }
            else {
                if (z_or_t) { //t
                    t_counter += -8;
                    indexList[currentIndex] = t + t_counter;
                }
                else { //z
                    z_counter += -8;
                    indexList[currentIndex] = z + z_counter;
                }
            }
            currentIndex++;
            z_or_t = !z_or_t;
        }
        z_counter += -4; t_counter += -4;
        if (z_or_t) {
            indexList[currentIndex] = t + t_counter; currentIndex++;
            indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 2nd ring
        }
        else {
            indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 2nd ring
            indexList[currentIndex] = t + t_counter; currentIndex++;
        }

        indexList[currentIndex] = indexList[currentIndex - 2];  //2nd Degenerate triangle
        t_counter = indexList[currentIndex - 2] - t;  z_or_t = !z_or_t;
        currentIndex++;
        for (unsigned int i = 0; i < corner_point_amount; i++) { //3rd ring [:-1] (8)
            if (false) {
                indexList[currentIndex] = indexList[currentIndex - 2];
                t_counter = indexList[currentIndex - 2] - t;
            }
            else {
                if (z_or_t) { //t
                    t_counter += 8;
                    indexList[currentIndex] = t + t_counter;
                }
                else { //z
                    z_counter += 8;
                    indexList[currentIndex] = z + z_counter;
                }
            }
            currentIndex++;
            z_or_t = !z_or_t;
        }
        z_counter += 4; t_counter += 4;
        if (z_or_t) {
            indexList[currentIndex] = t + t_counter; currentIndex++;
            indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 3rd ring (4)
        }
        else {
            indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 3rd ring (4)
            indexList[currentIndex] = t + t_counter; currentIndex++;
        }
        indexList[currentIndex] = indexList[currentIndex - 2];  //3rd Degenerate triangle
        t_counter = indexList[currentIndex - 2] - t;  z_or_t = !z_or_t;
        currentIndex++;

        z_counter += -6; t_counter += -6;
        indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 4th ring (-8)
        indexList[currentIndex] = t + t_counter; currentIndex++;

        for (unsigned int i = 0; i < corner_point_amount; i++) { //4th ring [:-1] (-8)
            if (false) {
                indexList[currentIndex] = indexList[currentIndex - 2];
                t_counter = indexList[currentIndex - 2] - t;
            }
            else {
                if (z_or_t) { //t
                    t_counter += -8;
                    indexList[currentIndex] = t + t_counter;
                }
                else { //z
                    z_counter += -8;
                    indexList[currentIndex] = z + z_counter;
                }
            }
            currentIndex++;
            z_or_t = !z_or_t;
        }
        if (!f) return currentIndex;
        if (next_ring_is_larger) {
            z_counter += 4;
            indexList[currentIndex] = z + z_counter;
            currentIndex++;
        }
        else {
            t_counter += 4;
            indexList[currentIndex] = t + t_counter;
            currentIndex++;
        }
        return currentIndex;
    }

    //middle is always t and the next ring is always smaller
    unsigned int fillRingMiddle(unsigned int t, unsigned int z, unsigned int t_total, unsigned int z_total,  unsigned int currentIndex) { //currentIndex is for the SPECIFIC ring not the currentIndex for the whole
        indexList[currentIndex] = z; indexList[currentIndex + 1] = t; //initial 2 points for base of the first triangle
        currentIndex += 2;
        unsigned int t_counter = 2; unsigned int z_counter = 4;
        unsigned int corner_point_amount = (2 * (((z_total / 4) - 1) - 2)) + 3; //amount of triangles between degens
        indexList[currentIndex] = z + z_counter; currentIndex++; //initial 2 points for 1st ring
        indexList[currentIndex] = t + t_counter; currentIndex++;

        bool z_or_t = 0; //0 = z 1 = t
        unsigned int degenCounter = 0; //start on zero since technically first triangle hasnt been drawn yet
        unsigned int degenCounterMax = 9999; //logic for setting this comes later
        for (unsigned int i = 0; i < corner_point_amount; i++) { //first ring [1:]
            if (z_or_t) { //t
                t_counter += 2;
                indexList[currentIndex] = t + t_counter;
            }
            else { //z
                z_counter += 8;
                indexList[currentIndex] = z + z_counter;
            }
            currentIndex++;
            z_or_t = !z_or_t;
        }



        indexList[currentIndex] = indexList[currentIndex - 2];  //First Degenerate triangle
        t_counter = indexList[currentIndex - 2] - t;  z_or_t = !z_or_t;
        currentIndex++;

        z_counter += -6; t_counter += -3;
        if (z_or_t) {
            indexList[currentIndex] = t + t_counter; currentIndex++;
            indexList[currentIndex] = z + z_counter; currentIndex++; //initial 2 points for 2nd ring (4)
        }
        else {
            indexList[currentIndex] = z + z_counter; currentIndex++; //initial 2 points for 2nd ring (4)
            indexList[currentIndex] = t + t_counter; currentIndex++;
        }

        for (unsigned int i = 0; i < corner_point_amount - 2; i++) { //2nd ring [1:-1] 
            if (z_or_t) { //t
                t_counter += -4;
                indexList[currentIndex] = t + t_counter;
            }
            else { //z
                z_counter += -8;
                indexList[currentIndex] = z + z_counter;
            }
            currentIndex++;
            z_or_t = !z_or_t;
        }

        z_counter += -4; t_counter += -2;
        if (z_or_t) {
            indexList[currentIndex] = t + t_counter; currentIndex++;
            indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 2nd ring
        }
        else {
            indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 2nd ring
            indexList[currentIndex] = t + t_counter; currentIndex++;
        }

        indexList[currentIndex] = indexList[currentIndex - 2];  //2nd Degenerate triangle
        t_counter = indexList[currentIndex - 2] - t;  z_or_t = !z_or_t;
        currentIndex++;
        for (unsigned int i = 0; i < corner_point_amount; i++) { //3rd ring [:-1] (8)
            if (z_or_t) { //t
                t_counter += 4;
                indexList[currentIndex] = t + t_counter;
            }
            else { //z
                z_counter += 8;
                indexList[currentIndex] = z + z_counter;
            }
            currentIndex++;
            z_or_t = !z_or_t;
        }

        z_counter += 4; t_counter += 2;
        if (z_or_t) {
            indexList[currentIndex] = t + t_counter; currentIndex++;
            indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 3rd ring (4)
        }
        else {
            indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 3rd ring (4)
            indexList[currentIndex] = t + t_counter; currentIndex++;
        }

        indexList[currentIndex] = indexList[currentIndex - 2];  //3rd Degenerate triangle
        t_counter = indexList[currentIndex - 2] - t;  z_or_t = !z_or_t;
        currentIndex++;

        z_counter += -6; t_counter += -3;
        indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 4th ring (-8)
        indexList[currentIndex] = t + t_counter; currentIndex++;


        for (unsigned int i = 0; i < corner_point_amount; i++) { //4th ring [:-1] (-8)
            if (z_or_t) { //t
                t_counter += -4;
                indexList[currentIndex] = t + t_counter;
            }
            else { //z
                z_counter += -8;
                indexList[currentIndex] = z + z_counter;
            }
            currentIndex++;
            z_or_t = !z_or_t;
        }

        if (true) {
            indexList[currentIndex] = z + z_counter + 2;
            currentIndex++;
            indexList[currentIndex] = z + z_counter;
            currentIndex++;
        }
        else {
            t_counter += 2;
            indexList[currentIndex] = t + t_counter;
            currentIndex++;
        }

        return currentIndex;
    }

    //middle is always z and the next ring is always smaller
    unsigned int fillRingToMiddle(unsigned int t, unsigned int z, unsigned int t_total, unsigned int z_total, unsigned int currentIndex) { //currentIndex is for the SPECIFIC ring not the currentIndex for the whole
        indexList[currentIndex] = z; indexList[currentIndex + 1] = t; //initial 2 points for base of the first triangle
        currentIndex += 2;
        unsigned int t_counter = 4; unsigned int z_counter = 2;
        unsigned int corner_point_amount = (2 * (((z_total / 4) - 1) - 2)) + 3; //amount of triangles between degens
        indexList[currentIndex] = z + z_counter; currentIndex++; //initial 2 points for 1st ring
        indexList[currentIndex] = t + t_counter; currentIndex++;

        bool z_or_t = 0; //0 = z 1 = t
        unsigned int degenCounter = 0; //start on zero since technically first triangle hasnt been drawn yet
        unsigned int degenCounterMax = 9999; //logic for setting this comes later
        for (unsigned int i = 0; i < corner_point_amount; i++) { //first ring [1:]
            if (z_or_t) { //t
                t_counter += 8;
                indexList[currentIndex] = t + t_counter;
            }
            else { //z
                z_counter += 4;
                indexList[currentIndex] = z + z_counter;
            }
            currentIndex++;
            z_or_t = !z_or_t;
        }


        indexList[currentIndex] = indexList[currentIndex - 2];  //First Degenerate triangle
        t_counter = indexList[currentIndex - 2] - t;  z_or_t = !z_or_t;
        currentIndex++;

        z_counter += -3; t_counter += -6;
        if (z_or_t) {
            indexList[currentIndex] = t + t_counter; currentIndex++;
            indexList[currentIndex] = z + z_counter; currentIndex++; //initial 2 points for 2nd ring (4)
        }
        else {
            indexList[currentIndex] = z + z_counter; currentIndex++; //initial 2 points for 2nd ring (4)
            indexList[currentIndex] = t + t_counter; currentIndex++;
        }

        for (unsigned int i = 0; i < corner_point_amount - 2; i++) { //2nd ring [1:-1] 
            if (z_or_t) { //t
                t_counter += -8;
                indexList[currentIndex] = t + t_counter;
            }
            else { //z
                z_counter += -4;
                indexList[currentIndex] = z + z_counter;
            }
            currentIndex++;
            z_or_t = !z_or_t;
        }

        z_counter += -2; t_counter += -4;
        if (z_or_t) {
            indexList[currentIndex] = t + t_counter; currentIndex++;
            indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 2nd ring
        }
        else {
            indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 2nd ring
            indexList[currentIndex] = t + t_counter; currentIndex++;
        }

        indexList[currentIndex] = indexList[currentIndex - 2];  //2nd Degenerate triangle
        t_counter = indexList[currentIndex - 2] - t;  z_or_t = !z_or_t;
        currentIndex++;
        for (unsigned int i = 0; i < corner_point_amount; i++) { //3rd ring [:-1] (8)
            if (z_or_t) { //t
                t_counter += 8;
                indexList[currentIndex] = t + t_counter;
            }
            else { //z
                z_counter += 4;
                indexList[currentIndex] = z + z_counter;
            }
            currentIndex++;
            z_or_t = !z_or_t;
        }

        z_counter += 2; t_counter += 4;
        if (z_or_t) {
            indexList[currentIndex] = t + t_counter; currentIndex++;
            indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 3rd ring (4)
        }
        else {
            indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 3rd ring (4)
            indexList[currentIndex] = t + t_counter; currentIndex++;
        }

        indexList[currentIndex] = indexList[currentIndex - 2];  //3rd Degenerate triangle
        t_counter = indexList[currentIndex - 2] - t;  z_or_t = !z_or_t;
        currentIndex++;

        z_counter += -3; t_counter += -6;
        indexList[currentIndex] = z + z_counter; currentIndex++; //last 2 points for 4th ring (-8)
        indexList[currentIndex] = t + t_counter; currentIndex++;


        for (unsigned int i = 0; i < corner_point_amount; i++) { //4th ring [:-1] (-8)
            if (z_or_t) { //t
                t_counter += -8;
                indexList[currentIndex] = t + t_counter;
            }
            else { //z
                z_counter += -4;
                indexList[currentIndex] = z + z_counter;
            }
            currentIndex++;
            z_or_t = !z_or_t;
        }

        if (true) {
            z_counter += 4;
            indexList[currentIndex] = z + z_counter;
            currentIndex++;
        }
        else {
            t_counter += 4;
            indexList[currentIndex] = t + t_counter;
            currentIndex++;
        }

        return currentIndex;
    }


    void sphereIndicesGeneratorV3() {
        unsigned int currentIndexAdjusted = currentIndex / VERTEX_SIZE;
        unsigned int t = 10 + currentIndexAdjusted; unsigned int t_total = 8; //current ring
        unsigned int z = 26 + currentIndexAdjusted; unsigned int z_total = 12; //what its connecting to
        unsigned int currentRingIndex = indiceIndex; //this is for each individual ring not for the entire shape
        unsigned int counter = 8 * 3;
        
        
        indexList[currentRingIndex] = 2 + currentIndexAdjusted; currentRingIndex++; //1->4
        indexList[currentRingIndex] = 8 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 0 + currentIndexAdjusted; currentRingIndex++; 
        
        indexList[currentRingIndex] = 4 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 6 + currentIndexAdjusted; currentRingIndex++;
        
        indexList[currentRingIndex] = 6 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 0 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 6 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 2 + currentIndexAdjusted; currentRingIndex++;
        
        indexList[currentRingIndex] = 10 + currentIndexAdjusted; currentRingIndex++; //assuming that theres a ring between [1](ring of 4 vertices) and middle
        indexList[currentRingIndex] = 2 + currentIndexAdjusted; currentRingIndex++; //4->8
        indexList[currentRingIndex] = 18 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 8 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 24 + currentIndexAdjusted; currentRingIndex++; //where it actually is
        indexList[currentRingIndex] = 8 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 20 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 4 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 12 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 4 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 16 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 6 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 22 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 6 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 14 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 2 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 10 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = 14 + currentIndexAdjusted; currentRingIndex++;
        
        currentRingIndex = fillRing(t, z, t_total, z_total, true, currentRingIndex, true);
        
        t_total += 4; z_total += 4;
        for (unsigned int i = 4; i < halfRingAmount; i++) { //approaching middle
            
            
            t = z;
            z += counter;
            
            counter += 8; //8 due to mirroring
            currentRingIndex = fillRing(t, z, t_total, z_total, true, currentRingIndex, true);
            t_total += 4; z_total += 4;
        }   
        //non middle-> middle
        t = z;
        z += counter;

        currentRingIndex = fillRingToMiddle(t, z, t_total, z_total, currentRingIndex);
        t_total += 4; z_total -= 4;
        //middle-> non middle


        t = z;
        z -= counter;
        z += 1; //to get to other side of mirror
        counter -= 8; //8 due to mirroring
        currentRingIndex = fillRingToMiddle(z, t, z_total, t_total, currentRingIndex);
        t_total -= 4; z_total -= 4;

        for (unsigned int i = halfRingAmount - 1; i > 2; i--) { //last ones
            t = z;
            z -= counter;
            counter -= 8; //8 due to mirroring
            currentRingIndex = fillRing(z, t, z_total, t_total, true, currentRingIndex, false);
            t_total -= 4; z_total -= 4;
        }

        t = z;
        z -= counter;
        counter -= 8; //8 due to mirroring
        //currentRingIndex = fillRing(t, z, t_total, z_total, false, currentRingIndex, false);
        //t_total -= 4; z_total -= 4;

        indexList[currentRingIndex] = t; currentRingIndex++; //assuming that theres a ring between [1](ring of 4 vertices) and middle
        indexList[currentRingIndex] = z; currentRingIndex++;
        indexList[currentRingIndex] = t+8; currentRingIndex++;
        indexList[currentRingIndex] = z+6; currentRingIndex++;
        indexList[currentRingIndex] = t+14; currentRingIndex++;
        indexList[currentRingIndex] = z+6; currentRingIndex++;
        indexList[currentRingIndex] = t+10; currentRingIndex++;
        indexList[currentRingIndex] = z+2; currentRingIndex++;
        indexList[currentRingIndex] = t+2; currentRingIndex++;
        indexList[currentRingIndex] = z+2; currentRingIndex++;
        indexList[currentRingIndex] = t+6; currentRingIndex++;
        indexList[currentRingIndex] = z+4; currentRingIndex++;
        indexList[currentRingIndex] = t+12; currentRingIndex++;
        indexList[currentRingIndex] = z+4; currentRingIndex++;
        indexList[currentRingIndex] = t+4; currentRingIndex++;
        indexList[currentRingIndex] = z; currentRingIndex++;
        indexList[currentRingIndex] = t; currentRingIndex++;
        indexList[currentRingIndex] = z+4; currentRingIndex++;
        
        indexList[currentRingIndex] = z; currentRingIndex++;
        indexList[currentRingIndex] = z+6; currentRingIndex++;
        indexList[currentRingIndex] = 1 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = z+2; currentRingIndex++;
        indexList[currentRingIndex] = z+4; currentRingIndex++;
        indexList[currentRingIndex] = z+4; currentRingIndex++;
        indexList[currentRingIndex] = 1 + currentIndexAdjusted; currentRingIndex++;
        indexList[currentRingIndex] = z+4; currentRingIndex++;
        indexList[currentRingIndex] = z; currentRingIndex++;
        std::cout << "CURRENTIndex: " << currentRingIndex << std::endl;

    }


    inline void drawPoint(unsigned int xIndex, unsigned int yIndex, unsigned int zIndex, float* holder, unsigned int &totalVertices, unsigned int startIndex) {
        holder[((0 + totalVertices) * VERTEX_SIZE) + 0 + startIndex] = uniqueValues[xIndex] + x;  //x 
        holder[((0 + totalVertices) * VERTEX_SIZE) + 1 + startIndex] = uniqueValues[yIndex] + y;  //y
        holder[((0 + totalVertices) * VERTEX_SIZE) + 2 + startIndex] = uniqueValues[zIndex] + z;  //z //point behind center
        //holder[((0 + totalVertices) * VERTEX_SIZE) + 3 + startIndex] = totalPoints;  //z //point behind center
        //holder[((0 + totalVertices) * VERTEX_SIZE) + 4 + startIndex] = totalPoints;  //z //point behind center
        //holder[((0 + totalVertices) * VERTEX_SIZE) + 5 + startIndex] = totalPoints;  //z //point behind center
        totalVertices++;
        totalPoints += 1.0 / (length / 7.0);
        //x,y,and z are position of sphere to offset by while the indexes are for indexing uniqueValues
    }

    inline void drawPointYMirrored(unsigned int xIndex, unsigned int yIndex1, unsigned int yIndex2, unsigned int zIndex, float* holder, unsigned int &totalVertices, unsigned int startIndex) {
        drawPoint(xIndex, yIndex1, zIndex, holder, totalVertices, startIndex);
        drawPoint(xIndex, yIndex2, zIndex, holder, totalVertices, startIndex);
    } 

    void sphereVerticesGeneratorV2(GLfloat* holder, unsigned int startIndex) {
        holder[(0 * VERTEX_SIZE) + startIndex] = x;
        holder[(0 * VERTEX_SIZE) + 1 + startIndex] = y + radius;
        holder[(0 * VERTEX_SIZE) + 2 + startIndex] = z;

        holder[(1 * VERTEX_SIZE) + startIndex] = x;
        holder[(1 * VERTEX_SIZE) + 1 + startIndex] = y - radius;
        holder[(1 * VERTEX_SIZE) + 2 + startIndex] = z;

        unsigned int totalVertices = 2;
        unsigned int ringAmountLess = ringAmount - 1;
        for (unsigned int ring = 1; ring < halfRingAmount; ring++) { //from [1] to [middle - 1] //everything is doubled to mirror across y
            //std::cout << ring << "            RING" << std::endl;
            drawPointYMirrored(halfRingAmount, ring, ringAmountLess - ring, halfRingAmount + ring, holder, totalVertices, startIndex); //x //point right of center //z //point behind center
            drawPointYMirrored(halfRingAmount, ring, ringAmountLess - ring, halfRingAmount - ring, holder, totalVertices, startIndex); //x //point right of center //z //point infront of center

            unsigned int extra = 1;
            for (; extra < ring; extra++) {
                drawPointYMirrored(halfRingAmount + extra, ring, ringAmountLess - ring, halfRingAmount + (ring - extra), holder, totalVertices, startIndex); //x //point right of center //z //point behind center
                drawPointYMirrored(halfRingAmount + extra, ring, ringAmountLess - ring, halfRingAmount - (ring - extra), holder, totalVertices, startIndex); //x //point right of center //z //point infront of center
                drawPointYMirrored(halfRingAmount - extra, ring, ringAmountLess - ring, halfRingAmount + (ring - extra), holder, totalVertices, startIndex); //x //point left of center //z //point behind center
                drawPointYMirrored(halfRingAmount - extra, ring, ringAmountLess - ring, halfRingAmount - (ring - extra), holder, totalVertices, startIndex); //x //point left of center //z //point infront of center
            } //2nd iteration of this loop for z 
                //extra logic of tips of the ring since they arent mirrored across the z
                drawPointYMirrored(halfRingAmount + extra, ring, ringAmountLess - ring, halfRingAmount + (ring - extra), holder, totalVertices, startIndex); //x //point rightmost of center
                drawPointYMirrored(halfRingAmount - extra, ring, ringAmountLess - ring, halfRingAmount - (ring - extra), holder, totalVertices, startIndex); //x //point leftmost of center
        }
        //extra logic is here for the middle ring since it doesnt mirror accross the z
        drawPoint(halfRingAmount, halfRingAmount, ringAmountLess, holder, totalVertices, startIndex);
        drawPoint(halfRingAmount, halfRingAmount, 0, holder, totalVertices, startIndex);
        
        unsigned int extra = 1;
        for (; extra < halfRingAmount; extra++) {
            drawPoint(halfRingAmount + extra, halfRingAmount, halfRingAmount + (halfRingAmount - extra), holder, totalVertices, startIndex); //x //point right of center //z //point behind center
            
            drawPoint(halfRingAmount + extra, halfRingAmount, halfRingAmount - (halfRingAmount - extra), holder, totalVertices, startIndex); //x //point right of center //z //point infront of center
            drawPoint(halfRingAmount - extra, halfRingAmount, halfRingAmount - (halfRingAmount - extra), holder, totalVertices, startIndex); //x //point left of center //z //point behind center
            drawPoint(halfRingAmount - extra, halfRingAmount, halfRingAmount + (halfRingAmount - extra), holder, totalVertices, startIndex); //x //point left of center //z //point infront of center
        }
        //extra logic of tips of the ring since they arent mirrored across the z
        drawPoint(halfRingAmount + extra, halfRingAmount, halfRingAmount, holder, totalVertices, startIndex); //x //point rightmost of center
        drawPoint(halfRingAmount - extra, halfRingAmount, halfRingAmount, holder, totalVertices, startIndex); //x //point leftmost of center
    }


    void uniqueValueGenerator(float* uniqueValues, unsigned int ringAmount, unsigned int halfRingAmount) { //uniqueValues is assumed to already be allocated

        uniqueValues[0] = radius;
        unsigned int indexer = 1;
        for (float i = 1; i < ringAmount; i++) {
            if (indexer == halfRingAmount) {
                uniqueValues[indexer] = 0.0;
            }
            else {
                uniqueValues[indexer] = radius * cos((M_PI * i) / (ringAmount - 1.0));
            }
            indexer++; //this is here because you can not index by a float
        }
    }

    float radius;
    //how many circles in a sphere
    unsigned int halfRingAmount;
    unsigned int ringAmount;
    unsigned int* verticeRingAmounts;
    //multiply these by radius of ring their applied to
    float* uniqueValues = nullptr;


    float totalPoints = 0;
    float incrementer = 0;

};

//class for making spheres, takes amount of rings as input to calculate vertices needed and etc
class SphereConstructor {
public:
    SphereConstructor(GLenum input_draw_type, unsigned int ringAmount, float X, float Y, float Z, float radius, Vector4 Color) {
        this->X = X;
        this->Y = Y;
        this->Z= Z;
        this->radius = radius;
        this->Color = Color;
        this->input_draw_type = input_draw_type;
        if (ringAmount < 8) { //makes sure there are atleast 5 rings
            ringAmount = 9;
        }
        else if (ringAmount % 2 == 0) { //makes sure its an odd amount of rings (goes down to prevent overflow and mins at 5)
            ringAmount -= 1;
        }
        this->ringAmount = ringAmount;
        halfRingAmount = ringAmount / 2; //this doesnt include the center ring
        verticeRingAmounts = new unsigned int[ringAmount]; //free after sphere is made
        calcRingAmounts(ringAmount, halfRingAmount, verticeRingAmounts); //fills in the vertice amount per ring (including 1st and last)
        length = calcLength(halfRingAmount, verticeRingAmounts);
        indexLength = calcIndexLength(halfRingAmount);
    }

    Sphere getSphere(void) {
        Sphere mySphere(input_draw_type, length, indexLength, ringAmount, verticeRingAmounts, halfRingAmount, radius, X, Y, Z, Color);
        return mySphere;
    }
    Sphere* mySphere = nullptr; //this is a ptr as to remain uninitialized
    unsigned int ringAmount;
    unsigned int halfRingAmount;
    unsigned int* verticeRingAmounts;
    unsigned int length;
    unsigned int indexLength;
    float X;
    float Y;
    float Z;
    float radius;
    Vector4 Color;
    GLenum input_draw_type;
protected:
    unsigned int calcLength(unsigned int halfRingAmount, unsigned int* verticeRingAmounts) {
        return static_cast<unsigned int>(((halfRingAmount * 0.5) * (verticeRingAmounts[halfRingAmount - 1])) * 2 + verticeRingAmounts[halfRingAmount] + 2); //fix this
    }
    unsigned int calcIndexLength(unsigned int halfRingAmount) {
        halfRingAmount -= 2;
        return 54 + halfRingAmount * (52 + (halfRingAmount - 1) * 8); //this is the doubled version
    }
    void calcRingAmounts(unsigned int ringAmount, unsigned int halfRingAmount, unsigned int* verticeRingAmounts) {
        verticeRingAmounts[0] = 1;
        verticeRingAmounts[1] = 4;
        verticeRingAmounts[ringAmount - 2] = 4;
        verticeRingAmounts[ringAmount - 1] = 1;
        verticeRingAmounts[halfRingAmount] = 4 * (halfRingAmount); //middle ring
        unsigned int i = 2;
        for (; i < halfRingAmount; i++) {
            int amt = 4 * i;
            verticeRingAmounts[i] = amt;
            verticeRingAmounts[ringAmount - (i + 1)] = amt;
        }
    }
};

//seperate class because the difference in drawing style demands more indices on the indexbuffer
class Filled_Rectangle3D : public Polygon3D {
public:
    Filled_Rectangle3D(GLenum input_draw_type, float X, float Y, float Z, float Height, float Width, float Depth, Vector4 Color = DEFAULT_COLOR) : Polygon3D(input_draw_type, 8, 36, X, Y, Z, Color) {
        height = Height;
        width = Width;
        depth = Depth;

        //called manually to prevent need for initialization check during runtime
        rectangleVerticesGenerator3D(parentList, currentIndex);
        //localList is a backup to preserve unrotated position (this preserves potential loss of data from floating point errors)
        copyToLocalList();
        //since order doesnt change we only need to generate indices once (if its needed in the future it will be made public for that usecase)
        rectangle3dIndexGenerator();
        localListInitialized = true;
        //drawRectangle();
    }

    //this function is assumed to run after localist has been initialized by copying from the parent list
    void drawRectangle() {
        //generate seperate list of vertices for each
        rectangleVerticesGenerator3D(parentList, currentIndex);
        //localList is a backup to preserve unrotated position (this preserves potential loss of data from floating point errors)
        rectangleVerticesGenerator3D(localList, 0);
        updateRotation();

    }

    void setXYZ(float X, float Y, float Z) {
        x = X;
        y = Y;
        z = Z;
        drawRectangle();
    }
    void setX(float X) {
        x = X;
        drawRectangle();
    }
    void setY(float Y) {
        y = Y;
        drawRectangle();
    }
    void setZ(float Z) {
        z = Z;
        drawRectangle();
    }

    //color is public to be modified the same way that the rotation is
    void update(void) {
        drawRectangle();
    }

    void setPositionAndColor(float X, float Y, float Z, Vector4 Color) {
        x = X;
        y = Y;
        z = Z;
        color = Color;
        drawRectangle();
    }

    //this version uses the indexBuffer 
    void drawSelfIndex(void) {
        if (drawType) {
            if (ready) {
                glDrawElements(GL_TRIANGLES, len, GL_UNSIGNED_INT, indexList + indiceIndex);
            }
        }
        else std::cout << "something odd happened line: " << __LINE__ << " of file: " << __FILE__ << std::endl;
    }
    void drawSelfIndex(bool iboBounded) {
        if (drawType) {
            if (!iboBounded) {
                glDrawElements(GL_TRIANGLES, len, GL_UNSIGNED_INT, indexList + indiceIndex);
            }
            else {
                glDrawElements(GL_TRIANGLES, len, GL_UNSIGNED_INT, (void*)(indiceIndex * sizeof(GLuint)));
            }
        }
        else std::cout << "something odd happened line: " << __LINE__ << " of file: " << __FILE__ << std::endl;
    }
    //chaningpoint is the currentindex of the element removed

protected:
    //holder is the destination
    //everything was hardcoded into it as constants for runtime efficency (its faster to read than generate new values)
    void rectangleVerticesGenerator3D(GLfloat* holder, unsigned int startIndex) {
        //x values
        //holder[0] = 1;
        holder[(0 * VERTEX_SIZE) + 0 + startIndex] = x + (width / 2);
        holder[(1 * VERTEX_SIZE) + 0 + startIndex] = x + (width / 2);
        holder[(2 * VERTEX_SIZE) + 0 + startIndex] = x - (width / 2);
        holder[(3 * VERTEX_SIZE) + 0 + startIndex] = x - (width / 2);
        holder[(4 * VERTEX_SIZE) + 0 + startIndex] = x + (width / 2);
        holder[(5 * VERTEX_SIZE) + 0 + startIndex] = x + (width / 2);
        holder[(6 * VERTEX_SIZE) + 0 + startIndex] = x - (width / 2);
        holder[(7 * VERTEX_SIZE) + 0 + startIndex] = x - (width / 2);
        //y values
        holder[(0 * VERTEX_SIZE) + 1 + startIndex] = y - (height / 2);
        holder[(1 * VERTEX_SIZE) + 1 + startIndex] = y + (height / 2);
        holder[(2 * VERTEX_SIZE) + 1 + startIndex] = y - (height / 2);
        holder[(3 * VERTEX_SIZE) + 1 + startIndex] = y + (height / 2);
        holder[(4 * VERTEX_SIZE) + 1 + startIndex] = y - (height / 2);
        holder[(5 * VERTEX_SIZE) + 1 + startIndex] = y + (height / 2);
        holder[(6 * VERTEX_SIZE) + 1 + startIndex] = y - (height / 2);
        holder[(7 * VERTEX_SIZE) + 1 + startIndex] = y + (height / 2);
        //z values
        holder[(0 * VERTEX_SIZE) + 2 + startIndex] = z - (depth / 2);
        holder[(1 * VERTEX_SIZE) + 2 + startIndex] = z - (depth / 2);
        holder[(2 * VERTEX_SIZE) + 2 + startIndex] = z - (depth / 2);
        holder[(3 * VERTEX_SIZE) + 2 + startIndex] = z - (depth / 2);
        holder[(4 * VERTEX_SIZE) + 2 + startIndex] = z + (depth / 2);
        holder[(5 * VERTEX_SIZE) + 2 + startIndex] = z + (depth / 2);
        holder[(6 * VERTEX_SIZE) + 2 + startIndex] = z + (depth / 2);
        holder[(7 * VERTEX_SIZE) + 2 + startIndex] = z + (depth / 2);
        //r values
        holder[(0 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(1 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(2 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(3 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(4 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(5 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(6 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        holder[(7 * VERTEX_SIZE) + 3 + startIndex] = color.x;
        //g values
        holder[(0 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(1 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(2 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(3 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(4 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(5 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(6 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        holder[(7 * VERTEX_SIZE) + 4 + startIndex] = color.y;
        //b values
        holder[(0 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(1 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(2 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(3 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(4 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(5 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(6 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        holder[(7 * VERTEX_SIZE) + 5 + startIndex] = color.z;
        //opacity values
        holder[(0 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(1 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(2 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(3 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(4 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(5 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(6 * VERTEX_SIZE) + 6 + startIndex] = color.w;
        holder[(7 * VERTEX_SIZE) + 6 + startIndex] = color.w;

    }
    void rectangle3dIndexGenerator(void) {
        //series of lines tells the GPU the order to connect points so we dont have to send duplicate vertices
        GLuint tempIndexList[36] = {
            5 + indiceIndex, 4 + indiceIndex, 6 + indiceIndex, //front
            6 + indiceIndex, 7 + indiceIndex, 5 + indiceIndex,
            5 + indiceIndex, 1 + indiceIndex, 4 + indiceIndex, //right
            4 + indiceIndex, 0 + indiceIndex, 1 + indiceIndex,
            1 + indiceIndex, 5 + indiceIndex, 7 + indiceIndex, //top
            7 + indiceIndex, 3 + indiceIndex, 1 + indiceIndex, 
            1 + indiceIndex, 3 + indiceIndex, 0 + indiceIndex, //back
            0 + indiceIndex, 3 + indiceIndex, 2 + indiceIndex,
            2 + indiceIndex, 6 + indiceIndex, 7 + indiceIndex, //left
            7 + indiceIndex, 3 + indiceIndex, 6 + indiceIndex,
            6 + indiceIndex, 4 + indiceIndex, 0 + indiceIndex, //bottom
            0 + indiceIndex, 2 + indiceIndex, 6 + indiceIndex,
        };
        memcpy(indexList + indiceIndex, tempIndexList, sizeof(GLuint) * 36);
    }
    float height;
    float width;
    float depth;
    bool filled;
};

struct windowData {
    int width;
    int height;
    int x;
    int y;
    int z;
    double mouseX;
    double mouseY;
    Quaternion camRotation;
    bool firstMouseInput = true;
    Vector3 camPosition = {0.0, 0.0, 0.0};
    std::vector<int> heldKeys;
};
//executes when the size of the window changes
void frame_buffer_size_callback(GLFWwindow* window, int width, int height) {
    windowData& windowDataObj = *static_cast<windowData*>(glfwGetWindowUserPointer(window));
    windowDataObj.width = width;
    windowDataObj.height = height;
    glViewport(0, 0, width, height);
    std::cout << width << " : width" << "\n" << height << " : height" << std::endl;
}

void mouse_input_callback(GLFWwindow* window, double x, double y) {
    windowData& windowDataObj = *static_cast<windowData*>(glfwGetWindowUserPointer(window));
    double deltaY = (y - windowDataObj.mouseY);
    double deltaX = (x - windowDataObj.mouseX);
    float screenWidth = static_cast<float>(windowDataObj.width);
    float screenHeight = static_cast<float>(windowDataObj.height);
    if (windowDataObj.firstMouseInput) {
        windowDataObj.firstMouseInput = false;
        return;
    }
    std::cout << deltaX << std::endl;
    if (deltaX != 0 && deltaY != 0 and false) {
        Quaternion rotator(cos((deltaX + deltaY) * 0.017 / 100.0) / 2, sin(deltaY * 0.017 / 100.0) / 2, sin(deltaX * 0.017 / 100.0) / 2, 0);
        rotator.normalize();
        windowDataObj.camRotation = windowDataObj.camRotation * rotator;
    } if (deltaX != 0) {
        Quaternion rotator(cos(deltaX * 0.017 / 100.0) / 2, 0, sin(deltaX * 0.017 / 100.0) / 2, 0);
        rotator.normalize();
        windowDataObj.camRotation = windowDataObj.camRotation * rotator;
    } if (deltaY != 0) {
        Quaternion rotator(cos(deltaY * 0.017 / 100.0) / 2, sin(deltaY * 0.017 / 100.0) / 2, 0, 0);
        rotator.normalize();
        windowDataObj.camRotation = rotator * windowDataObj.camRotation;
    }
    windowDataObj.mouseX = x;
    windowDataObj.mouseY = y;
}
/*
[in]	window	The window that received the event.
[in]	key	The keyboard key that was pressed or released.
[in]	scancode	The system-specific scancode of the key.
[in]	action	GLFW_PRESS, GLFW_RELEASE or GLFW_REPEAT. Future releases may add more actions.
[in]	mods	Bit field describing which modifier keys were held down.    */



void held_key_callback(GLFWwindow* window, int key) { //this function is for individual held keys
    switch (key) {
    case GLFW_KEY_A: {  //code canRotation.rotate(vector3)
        windowData& windowDataObj = *static_cast<windowData*>(glfwGetWindowUserPointer(window));
        Vector3 movementVector = windowDataObj.camRotation.getRight();
        movementVector = movementVector * -0.1;
        windowDataObj.camPosition.x += movementVector.x;
        windowDataObj.camPosition.y += movementVector.y;
        windowDataObj.camPosition.z += movementVector.z;
        break;
    } case GLFW_KEY_D: {
        windowData& windowDataObj = *static_cast<windowData*>(glfwGetWindowUserPointer(window));
        Vector3 movementVector = windowDataObj.camRotation.getRight();
        movementVector = movementVector * 0.1;
        windowDataObj.camPosition.x += movementVector.x;
        windowDataObj.camPosition.y += movementVector.y;
        windowDataObj.camPosition.z += movementVector.z;
        break;
    } case GLFW_KEY_S: {
        windowData& windowDataObj = *static_cast<windowData*>(glfwGetWindowUserPointer(window));
        Vector3 movementVector = windowDataObj.camRotation.getForward();
        movementVector = movementVector * -0.1;
        windowDataObj.camPosition.x += movementVector.x;
        windowDataObj.camPosition.y += movementVector.y;
        windowDataObj.camPosition.z += movementVector.z;
        break;
    } case GLFW_KEY_W: {
        windowData& windowDataObj = *static_cast<windowData*>(glfwGetWindowUserPointer(window));
        Vector3 movementVector = windowDataObj.camRotation.getForward();
        movementVector = movementVector * 0.1;
        windowDataObj.camPosition.x += movementVector.x;
        windowDataObj.camPosition.y += movementVector.y;
        windowDataObj.camPosition.z += movementVector.z;
        break;
    } case GLFW_KEY_E: {
        windowData& windowDataObj = *static_cast<windowData*>(glfwGetWindowUserPointer(window));
        Vector3 movementVector = windowDataObj.camRotation.getUp();
        movementVector = movementVector * 0.1;
        windowDataObj.camPosition.x += movementVector.x;
        windowDataObj.camPosition.y += movementVector.y;
        windowDataObj.camPosition.z += movementVector.z;
        break;
    } case GLFW_KEY_Q: {
        windowData& windowDataObj = *static_cast<windowData*>(glfwGetWindowUserPointer(window));
        Vector3 movementVector = windowDataObj.camRotation.getUp();
        movementVector = movementVector * -0.1;
        windowDataObj.camPosition.x += movementVector.x;
        windowDataObj.camPosition.y += movementVector.y;
        windowDataObj.camPosition.z += movementVector.z;
        break;
    }
    }
}

void held_keys_callback(GLFWwindow* window) { //this function loops through the held key list
    windowData& windowDataObj = *static_cast<windowData*>(glfwGetWindowUserPointer(window));
    for (unsigned int i = 0; i < windowDataObj.heldKeys.size(); i++) {
        int key = windowDataObj.heldKeys[i];
        held_key_callback(window, key);
    }
}


void keyboard_input_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    std::cout << "KEYBOARD INPUT CALLBACK";
    windowData& windowDataObj = *static_cast<windowData*>(glfwGetWindowUserPointer(window));
    auto iterator = std::find(windowDataObj.heldKeys.begin(), windowDataObj.heldKeys.end(), key); //the iterator is an "iterator" type
    if (iterator != windowDataObj.heldKeys.end()) {
        if (action == GLFW_RELEASE) windowDataObj.heldKeys.erase(iterator);
    } else {
        if (action == GLFW_PRESS) windowDataObj.heldKeys.push_back(key);
    }
    std::cout << std::endl;
}


int main(void)
{
    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    glfwSetFramebufferSizeCallback(window, frame_buffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_input_callback);
    glfwSetKeyCallback(window, keyboard_input_callback);
    //https://www.glfw.org/docs/3.0/group__input.html#ga7dad39486f2c7591af7fb25134a2501d documentation for mouse input callback
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    //Rectangle rect2(GL_LINE_LOOP, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f);
    //Rectangle3D rect4(GL_LINE_LOOP, 0.0f, 0.0f, 0.1f, 0.5f, 0.5f, 0.5f, { 1.0f, 1.0f, 1.0f, 1.0f });
    Rectangle3D rect3(GL_LINE_LOOP, 0.0f, 0.0f, -3.0f, 1.0f, 1.0f, 1.0f);
    Rectangle3D rect4(GL_LINE_LOOP, 0.0f, 0.0f, -5.0f, 1.0f, 1.0f, 1.0f);

    SphereConstructor constructor(GL_LINES, 24, 0.0f, 0.0f, -10.0f, 0.25f, { 1.0f, 1.0f, 1.0f, 1.0f });
    Sphere sphere3 = constructor.getSphere();
    
    //errors since its pre rect2


    //Rectangle rect(GL_LINE_LOOP, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f);
    //Rectangle rect3(GL_LINE_LOOP, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f);

    //Circle circle(GL_LINE_LOOP, 120, 0.4f, 0.5f, 0.5f, false);

    //rect2.~Rectangle();
    //new is added because stack based objects are inevitably removed so this prevents bugs from removing twice
    //new (&rect2) Rectangle(GL_LINE_LOOP, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f);

    GLfloat* vertices = Shape::getParentList();
    GLuint* indices = Shape::getIndexList();
    unsigned int verticeAmt = Shape::getParentListSize();
    unsigned int indexAmt = Shape::getIndexListSize();
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    }

    //glCompileShader turns the shader into readable machinecode for the gpu
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
    }


    GLuint shaderProgram = glCreateProgram();

    

    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);

    glLinkProgram(shaderProgram);
    //shaders are alr in the refrence object can be deleted
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    //cpu >> gpu is slow so stuff is sent in big batches(buffers)
    GLuint VAO, VBO, IBO;

    //normally an array of refrences
    glGenVertexArrays(1, &VAO);
    //VAO comes first.
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * verticeAmt, vertices, GL_STATIC_DRAW);
    Shape::setReady(true);

    //2nd int is amnt of values(vertices)
    //GL_FALSE is since our coords arent ints

    //starting index, amnt of points, dataType, normalize(may implement later), total size per point(VERTEX_SIZE), void ptr is just an offset(from start of VBO) for pointer arithematic
    //Position
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, VERTEX_SIZE * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    //Color (4th value is opacity)
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, VERTEX_SIZE * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    //total VERTEX_SIZE is summation of 2nd values for all of these


    glGenBuffers(1, &IBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * indexAmt, indices, GL_LINES);

    //glBindBuffer(GL_ARRAY_BUFFER, 0); //unbind vbo by binding to 0
    //glBindVertexArray(0); //unbind vao by binding to 0
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    /* Loop until the user closes the window */
    double frames = 0;
    float fps = 60;
    //rect.changeRotation(0.1f);
    float lastWidth, lastHeight;
    windowData windowDataObj;
    glfwSetWindowUserPointer(window, &windowDataObj);

    Quaternion rotator(-0.801143615547, 0.199490714701, 0.199490714701, 0.199490714701);
    rotator.normalize();
    double lastTime = glfwGetTime();
    double currentTime = glfwGetTime();
    SLERP newSlerp(windowDataObj.camRotation, rotator);
    double lastX = 0.0;
    double lastY = 0.0;
    double time = 0.0;

    


    GLint perspective_matrix_location = glad_glGetUniformLocation(shaderProgram, "u_perspectiveMatrix");
    GLint view_matrix_location = glad_glGetUniformLocation(shaderProgram, "u_viewMatrix");
    GLint position_matrix_location = glad_glGetUniformLocation(shaderProgram, "u_positionMatrix");

    while (!glfwWindowShouldClose(window))
    {   
        float normalizedSin = sin(static_cast<float>(glfwGetTime()));
        float normalizedCos = cos(static_cast<float>(glfwGetTime()));

        currentTime = glfwGetTime();
        //background color
        glClearColor(0.07f, 0.13f, 0.17f, 1.0f);
        float absSin = (normalizedSin / 2) + 0.5f;
        float absCos = (normalizedCos / 2) + 0.5f;
        //glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        //render here
        glClear(GL_COLOR_BUFFER_BIT);
        if (lastX != windowDataObj.mouseX || lastY != windowDataObj.mouseY) {
            //logic
            lastX = windowDataObj.mouseX;
            lastY = windowDataObj.mouseY;
            //rect3.setX(rect3.Z() * rect3.Z() * (((lastX / windowDataObj.width) - 0.5)));
            //rect3.setY(rect3.Z() * rect3.Z() * (-((lastY / windowDataObj.height) - 0.5)));
        }
        if ((currentTime - lastTime) > FRAME_RATE_DIVISOR) {
            //The fps loop (locked to 60fps)
            held_keys_callback(window);
            Quaternion oldQ1 = windowDataObj.camRotation;
            lastTime = currentTime;
            rotator.normalize();
            newSlerp.addTime(0.01);
            //rect3.color = { absCos , absSin , (absCos + absSin) / 2, rect3.color.w };
            //windowDataObj.camRotation = newSlerp.getQ3();
            //rect3.rotationQuaternion.normalize();
            rect3.updateRotation();
            rect3.update();
        }
        //rect.setXY(normalizedCos, normalizedSin);
        //circle.setRadius(normalizedSin);
        //vertices[0] = sin(currentFrame);
        //vertices[3] = sin(currentFrame);
        //vertices[6] = sin(currentFrame);
        //vertices[9] = sin(currentFrame);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(GLfloat) * verticeAmt, vertices);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * indexAmt, indices, GL_STATIC_DRAW);
        //Shape::setReady(true);

        glUseProgram(shaderProgram);

        float fov = 1.134464f;
        float f = 1.0 / tan(fov / 2.0f); //"   uniform mat4 uProjection;\n" //create the perspective matrix on the cpu then use a uniform to import to vertex shader
        float screenHeight = 1.0f;
        float screenWidth = 0.0f;
        float far_plane = 100.0f;
        float near_plane = 0.1f;

        windowData& windowDataObj = *static_cast<windowData*>(glfwGetWindowUserPointer(window));
        matrix3x3 rotation_matrix = windowDataObj.camRotation.getMatrix();
        rotation_matrix.m[1][1];
        

        float camera_matrix2[16] = { 
            rotation_matrix.m[0][0], rotation_matrix.m[1][0], rotation_matrix.m[2][0], 0,
            rotation_matrix.m[0][1], rotation_matrix.m[1][1], rotation_matrix.m[2][1], 0,
            rotation_matrix.m[0][2], rotation_matrix.m[1][2], rotation_matrix.m[2][2], 0,
            0, 0, 0, 1
        };


        float camera_matrix1[16] = { //combine 1 and 2
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            -windowDataObj.camPosition.x, -windowDataObj.camPosition.y, -windowDataObj.camPosition.z, 1
        };



        float camera_matrix3[16] = { //combine 1 and 2
            rotation_matrix.m[0][0], rotation_matrix.m[1][0], rotation_matrix.m[2][0], 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            rotation_matrix.m[0][0] * (-windowDataObj.camPosition.x) + rotation_matrix.m[0][1] * (-windowDataObj.camPosition.y) + rotation_matrix.m[0][2] * (-windowDataObj.camPosition.z), rotation_matrix.m[1][0] * (-windowDataObj.camPosition.x) + rotation_matrix.m[1][1] * (-windowDataObj.camPosition.y) + rotation_matrix.m[1][2] * (-windowDataObj.camPosition.z), rotation_matrix.m[2][0] * (-windowDataObj.camPosition.x) + rotation_matrix.m[2][1] * (-windowDataObj.camPosition.y) + rotation_matrix.m[2][2] * (-windowDataObj.camPosition.z), 0
        };



        float perspective_matrix[16] = {
            f / (static_cast<float>(windowDataObj.width) / static_cast<float>(windowDataObj.height)), 0, 0, 0,   // column 0
            0, f, 0, 0,   // column 1
            0, 0, -(far_plane + near_plane) / (far_plane - near_plane), -1.0f,   // column 2
            0, 0,  -(2.0f * far_plane * near_plane) / (far_plane - near_plane), 0    // column 3
        };

        glUniformMatrix4fv(perspective_matrix_location, 1, GL_FALSE, perspective_matrix);
        glUniformMatrix4fv(view_matrix_location, 1, GL_FALSE, camera_matrix2);
        glUniformMatrix4fv(position_matrix_location, 1, GL_FALSE, camera_matrix1);
        glBindVertexArray(VAO);
        //circle.drawSelf();
        //rect.changeRotation(deltaTime * normalizedSin * 10);
        //rect.drawSelf();
        sphere3.drawSelfIndex(true);
        //sphere3.drawSelf();
        //rect2.changeRotation(currentTime * normalizedSin * -10);
        //rect2.drawSelf();
        rect4.drawSelfIndex(true);
        //std::cout << "Q1w: " << rect3.rotationQuaternion.w << "   Q1x: " << rect3.rotationQuaternion.x << "    Q1y: " << rect3.rotationQuaternion.y << " Q1x: " << rect3.rotationQuaternion.z << std::endl;
        rect3.rotationQuaternion.normalize();
        rect3.drawSelfIndex(true);
        //rect3.drawSelf();
        //rect4.drawSelfIndex(true);

        //circle2.drawSelf();
        /* Swap front(output) and back(input) buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();
    }
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    //clearing

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glDeleteBuffers(1, &IBO);

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shaderProgram);

    glfwTerminate();
    return 0;
}