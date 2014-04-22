#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "./common.h"
#include "./bb.h"
#include "./mesh.h"
#include "./io.h"
#include "./material.h"
#include "./texture.h"

using namespace std;

// Added global Varibles, Refactor when possible

float lightG[3] = {100, 100, 100};
float lightGstar[3] = {100, 100, 100};
bool disp_scene = false;
bool disp_stargate = false;

float Scene1_Eye[3] = {1.66005e-08, 0.139173, 0.990268};
float Scene1_Up[3] = {-2.33307e-09, 0.990268, -0.139173};
float Scene1_Center[3] = {-0.199047, 9.67938, 23.8995};
float Scene1_EyeDis = 27.3469;

float eyeDist;    // Eye Distance from center
Vec3f eyeDir;     // Eye direction relative to Center
Vec3f center;     // Position in global space that the eye looks at
Vec3f upNorm;     // Vector that points up
float white[3] = {1.0f, 1.0f, 1.0f};
float black[3] = {0.0f, 0.0f, 0.0f};

static GLfloat lightPosition[4];
static GLfloat floorPlane[4] = {0.0f, 1.0f, 0.0f, 0.0f};
static GLfloat floorShadow[4][4];

float miny = 0;

Mesh mesh;

GLuint* texture_ids;

// window parameters
int window_width = 800, window_height = 600;
float window_aspect = window_width / static_cast<float>(window_height);

bool scene_lighting;

// takes in an eye location and converts to rotation and length
// Warning may only work with positive X Y Z


double sin_deg(double theta) {
  return sin(theta * M_PI / 180.0);
}
double cos_deg(double theta) {
  return cos(theta * M_PI / 180.0);
}


// rotates P around a by theta. P and axis must be normalized
Vec3f rotateArbitrary(float theta, Vec3f p, Vec3f a) {
  a = a.unit();
  float d = sqrt(a[1]*a[1] + a[2]*a[2]);
  // Rx( thetax )
  Vec3f rxx = Vec3f::makeVec(1, 0, 0);
  Vec3f rxy = Vec3f::makeVec(0, a[2]/d, -a[1]/d);
  Vec3f rxz = Vec3f::makeVec(0, a[1]/d, a[2]/d);
  // Ry( thetay )
  Vec3f ryx = Vec3f::makeVec(d, 0, -a[0]);
  Vec3f ryy = Vec3f::makeVec(0, 1, 0);
  Vec3f ryz = Vec3f::makeVec(a[0], 0, d);
  // Rz( theta )
  Vec3f rzx = Vec3f::makeVec(cos_deg(theta), -sin_deg(theta), 0);
  Vec3f rzy = Vec3f::makeVec(sin_deg(theta), cos_deg(theta), 0);
  Vec3f rzz = Vec3f::makeVec(0, 0, 1);
  p = Vec3f::makeVec(p*rxx, p*rxy, p*rxz);
  p = Vec3f::makeVec(p*ryx, p*ryy, p*ryz);
  p = Vec3f::makeVec(p*rzx, p*rzy, p*rzz);
  // Rx( -thetax )
  rxy[2] = -rxy[2];
  rxz[1] = -rxz[1];
  // Ry( -thetay )
  ryx[2] = -ryx[2];
  ryz[0] = -ryz[0];
  p = Vec3f::makeVec(p*ryx, p*ryy, p*ryz);
  p = Vec3f::makeVec(p*rxx, p*rxy, p*rxz);
  return p;
}

// Rotate Eye in degrees
// eyeLength must be initialized by initEye()
void rotateEye(double x, double y) {
  Vec3f eye = eyeDir.unit();
  Vec3f up = upNorm.unit();
  Vec3f right = (eye^up).unit();
  eye = rotateArbitrary(y, eye, right);
  Vec3f up2 = rotateArbitrary(y, up, right);
  eye = eye.unit();
  eye = rotateArbitrary(x, eye, up);
  eyeDir = eye.unit();
  upNorm = up2.unit();
  /*
  cout << "Eye:    " << eyeDir << endl;
  cout << "UP:     " << upNorm << endl;
  cout << "Center: " << center << endl;
  cout << "EyeDis: " << eyeDist << endl;
  */
}

void initLookup() {
  if (disp_scene) {
    eyeDir  = (Vec3f::makeVec(Scene1_Eye)).unit();
    center  = Vec3f::makeVec(Scene1_Center);
    upNorm  = Vec3f::makeVec(Scene1_Up).unit();
    eyeDist = Scene1_EyeDis;
  } else {
    eyeDir  = (Vec3f::makeVec(1, 0, 0)).unit();
    center  = Vec3f::makeVec(0, 0, 0);
    upNorm  = Vec3f::makeVec(0, 1, 0);
    eyeDist = 5.6;
    rotateEye(-45, 0);
  }
}

void adjustView() {
  Vec3f c = mesh.bb().center();
  Vec3f d = Vec3f::makeVec(mesh.bb().xdim(),
                           mesh.bb().ydim(),
                           mesh.bb().zdim());
  int length = d.norm()*1.4f;
  center = c;
  if (length > 0) eyeDist = length;
  miny = mesh.bb().min[1];
}


/*
   This code was borrowed from an internet guide called
   Shadows, Reflections, Lighting, Textures. Easy with OpenGL!
   It is largely unmodified and we make no claim to ownership
   http://www.opengl.org/archives/resources/code/samples/mjktips/TexShadowReflectLight.html
*/
enum {
  X, Y, Z, W
};
/* Create a matrix that will project the desired shadow. */
void shadowMatrix(GLfloat shadowMat[4][4],
  GLfloat groundplane[4],
  GLfloat lightpos[4]) {
  GLfloat dot;
  /* Find dot product between light position vector and ground plane normal. */
  dot = groundplane[X] * lightpos[X] +
    groundplane[Y] * lightpos[Y] +
    groundplane[Z] * lightpos[Z] +
    groundplane[W] * lightpos[W];

  shadowMat[0][0] = dot - lightpos[X] * groundplane[X];
  shadowMat[1][0] = 0.f - lightpos[X] * groundplane[Y];
  shadowMat[2][0] = 0.f - lightpos[X] * groundplane[Z];
  shadowMat[3][0] = 0.f - lightpos[X] * groundplane[W];

  shadowMat[X][1] = 0.f - lightpos[Y] * groundplane[X];
  shadowMat[1][1] = dot - lightpos[Y] * groundplane[Y];
  shadowMat[2][1] = 0.f - lightpos[Y] * groundplane[Z];
  shadowMat[3][1] = 0.f - lightpos[Y] * groundplane[W];

  shadowMat[X][2] = 0.f - lightpos[Z] * groundplane[X];
  shadowMat[1][2] = 0.f - lightpos[Z] * groundplane[Y];
  shadowMat[2][2] = dot - lightpos[Z] * groundplane[Z];
  shadowMat[3][2] = 0.f - lightpos[Z] * groundplane[W];

  shadowMat[X][3] = 0.f - lightpos[W] * groundplane[X];
  shadowMat[1][3] = 0.f - lightpos[W] * groundplane[Y];
  shadowMat[2][3] = 0.f - lightpos[W] * groundplane[Z];
  shadowMat[3][3] = dot - lightpos[W] * groundplane[W];
}
  //  Setups a light for opengl
void setupLight(GLenum light, Vec3f LightPos, Vec3f at) {
  Vec3f lightDir = (at-LightPos).unit();
  float pos[4] = {LightPos[0], LightPos[1], LightPos[2], 0};
  float dir[4] = {lightDir[0], lightDir[1], lightDir[2], 0};
  glLightfv(light, GL_POSITION, pos);
  glLightfv(light, GL_SPOT_DIRECTION, dir);
  glLightf(light, GL_SPOT_CUTOFF, 90.0);
  glLightfv(light, GL_AMBIENT, black);
  glLightfv(light, GL_DIFFUSE, white);
  glLightfv(light, GL_SPECULAR, white);
}

  // Renders a mesh
void render(const Mesh &m, bool text) {
  if (text) glEnable(GL_TEXTURE_2D);
  glCullFace(GL_BACK);
  for (int j = 0; j < m.polygons.size(); j++) {
    if (text) {
      int matIdx = m.polygon2material(j);
      Material mat = m.material(matIdx);
      glBindTexture(GL_TEXTURE_2D, mat.texture_id());
      glMaterialfv(GL_FRONT, GL_AMBIENT, mat.ambient().x);
      glMaterialfv(GL_FRONT, GL_DIFFUSE, mat.diffuse().x);
      glMaterialfv(GL_FRONT, GL_SPECULAR, mat.specular().x);
      glMaterialf(GL_FRONT, GL_SHININESS, mat.specular_coeff());
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_SRC_ALPHA);
    }
    glBegin(GL_POLYGON);
    for (int k = 0; k < m.polygons[j].size(); k++) {
      int vID = m.polygons[j][k];
      Vec3f v = m.vertices[vID].v;
      Vec3f n = m.vertices[vID].n;
      int tID = m.textPolygons[j][k];
      Vec3f t = m.texVertices[tID];
      glNormal3fv(n.x);
      if (text) glTexCoord2f(t[0], 1-t[1]);
      glVertex3fv(v.x);
    }
    glEnd();
  }
  glDisable(GL_CULL_FACE);
  if (text) glDisable(GL_TEXTURE_2D);
}

void DrawFloor(float W, float H, float w, float h) {
  float a = H/h, b = W/w;
  int M = static_cast<int>(floor(a+0.5f));
  int N = static_cast<int>(floor(b+0.5f));
  int i = 0, j = 0;
  Vec3f u = {w, 0, 0}, v = {0, 0, h}, r = {-(N/2)*w, 0, -(M/2)*h};
  Vec3f p0, p1, p2, p3;
  glEnable(GL_POLYGON_SMOOTH);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
  glBegin(GL_QUADS);
  glColor3f(0.7, 0.7, 0.7);
  for (j = 0; j < N; j++) {
    p0 = r+u*static_cast<float>(j%2);
    for (i = j%2; i < M; i += 2) {
      p1 = p0+u;
      p2 = p1+v;
      p3 = p2-u;

      glVertex3fv(p3.x);
      glVertex3fv(p2.x);
      glVertex3fv(p1.x);
      glVertex3fv(p0.x);

      p0 += u*2.0f;
    }
    r += v;
  }
  glEnd();
}

Vec3f scaleCenter(Vec3f v, float radius) {
  float a = sqrt(v[0]*v[0]+v[1]*v[1])/radius;
  v[2] *= 20*(1-a)*(1-a);
  return v;
}

Vec3f portalVec(float r, float t, float phi, float radius) {
  float x, y, z;
  x = r*cos_deg(t);
  y = r*sin_deg(t);
  float xx = x*x;
  float yy = y*y;
  z = .01*cos_deg(phi+2*360*sqrt(xx+yy)/radius);
  Vec3f c = Vec3f::makeVec(x, y, z);
  return scaleCenter(c, radius);
}

Vec3f portalNormal(float r, float t, float phi, float radius) {
  float rstep = radius/64;
  float tstep = 360/72.0f;
  Vec3f c = portalVec(r, t, phi, radius);
  Vec3f u = portalVec(r, t-tstep, phi, radius);
  Vec3f e = portalVec(r+rstep, t, phi, radius);
  Vec3f w = portalVec(r, t+tstep, phi, radius);
  Vec3f d = portalVec(r-rstep, t, phi, radius);
  /*
  c[2] /= .02;
  u[2] /= .02;
  e[2] /= .02;
  w[2] /= .02;
  d[2] /= .02;
  */
  Vec3f a1 = c - u;
  Vec3f b1 = c - e;
  Vec3f a2 = c - w;
  Vec3f b2 = c - d;
  Vec3f n1 = (b1.unit()^a1.unit()).unit();
  Vec3f n2 = (b2.unit()^a2.unit()).unit();
  Vec3f n = (n1+n2)/2;
  return Vec3f::zero()-n.unit();
}


void DrawPortal(float radius, float phi) {
  float rstep = radius/64;
  float tstep = 360/72;
  glColor4f(0.7, 0.9, 1.0, 1.0);
  float a[3] = {.1, .1, .1};
  float d[3] = {0.2, 0.2, 0.4};
  float s[3] = {0.6, 0.6, 1.0};
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, a);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, d);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, s);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30);
  for (float r = 0; r < radius; r += rstep) {
    for (float t = 0; t < 360; t += tstep) {
      glBegin(GL_POLYGON);
      Vec3f v;
      Vec3f n;
      v = portalVec(r, t+tstep, phi, radius);
      n = portalNormal(r, t+tstep, phi, radius);
      // v = scaleCenter(v, radius);
      glNormal3fv(n.x);
      glVertex3fv(v.x);
      v = portalVec(r+rstep, t+tstep, phi, radius);
      n = portalNormal(r+rstep, t+tstep, phi, radius);
      // v = scaleCenter(v, radius);
      glNormal3fv(n.x);
      glVertex3fv(v.x);
      v = portalVec(r+rstep, t, phi, radius);
      n = portalNormal(r+rstep, t, phi, radius);
      // v = scaleCenter(v, radius);
      glNormal3fv(n.x);
      glVertex3fv(v.x);
      v = portalVec(r, t, phi, radius);
      n = portalNormal(r, t, phi, radius);
      // v = scaleCenter(v, radius);
      glNormal3fv(n.x);
      glVertex3fv(v.x);
      glEnd();
    }
  }
}


void Display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // Set up View
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40.0, window_aspect, 1, 1500);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  Vec3f eye = eyeDir*eyeDist + center;
  gluLookAt(eye[0], eye[1], eye[2],
            center[0], center[1], center[2],
            upNorm[0], upNorm[1], upNorm[2]);
  Vec3f LightPos = eye;
  if (scene_lighting) LightPos = Vec3f::makeVec(lightG);
  if (disp_stargate) LightPos = Vec3f::makeVec(lightGstar);
  setupLight(GL_LIGHT0, LightPos, center);
  glEnable(GL_CULL_FACE);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_RESCALE_NORMAL);
  glPushMatrix();
  if (disp_stargate) {
    glRotatef(-90, 1, 0, 0);
  }
  render(mesh, true);
  glPopMatrix();


  // Cool looking portal
  if (disp_stargate) {
    glPushMatrix();
    static float p = 0;
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glTranslated(.05, -.1, 0);
    DrawPortal(.95, p);
    p += 20;
    glPopMatrix();
  }


  lightPosition[0] = LightPos[0];
  lightPosition[1] = LightPos[1];
  lightPosition[2] = LightPos[2];
  lightPosition[3] = 0;


  if (scene_lighting) {
    glPushMatrix();
    shadowMatrix(floorShadow, floorPlane, lightPosition);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glColor4f(0.0, 0.0, 0.0, 0.0);
    glMultMatrixf(reinterpret_cast<GLfloat *>(floorShadow));
    render(mesh, false);
    glDisable(GL_BLEND);
    glPopMatrix();
  }

  /* 
  glPushMatrix();
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
  glDisable(GL_BLEND);
  glTranslatef(0.0f, -0.04f, 0.0f);
  DrawFloor(800, 800, 80, 80);
  glPopMatrix();
  */
  glFlush();
  glutSwapBuffers();
  if (disp_stargate) {
    usleep(32000);
    glutPostRedisplay();
  }
}

void PrintMatrix(GLfloat* m) {
  cout.precision(2);
  int w = 6;
  for (int i = 0; i < 4; ++i) {
    cout << setprecision(2) << setw(w) << m[i] << " "
         << setprecision(2) << setw(w) << m[i+4] << " "
         << setprecision(2) << setw(w) << m[i+8] << " "
         << setprecision(2) << setw(w) << m[i+12] << " "
         << endl;
  }
  cout << endl;
}

void PrintMatrix(GLint matrix) {
  GLfloat m[16];
  glGetFloatv(matrix, m);
  PrintMatrix(m);
}

void PrintMatrix() {
  PrintMatrix(GL_MODELVIEW_MATRIX);
}

void LoadMatrix(GLfloat* m) {
  // transpose to column-major
  for (int i = 0; i < 4; ++i) {
    for (int j = i; j < 4; ++j) {
      swap(m[i*4+j], m[j*4+i]);
    }
  }
  glLoadMatrixf(m);
}

void MultMatrix(GLfloat* m) {
  // transpose to column-major
  for (int i = 0; i < 4; ++i) {
    for (int j = i; j < 4; ++j) {
      swap(m[i*4+j], m[j*4+i]);
    }
  }
  glMultMatrixf(m);
}

void Init() {
  initLookup();
  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);
  glDepthFunc(GL_LEQUAL);
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

  // resize the window
  window_aspect = window_width/static_cast<float>(window_height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(40.0, window_aspect, 1, 1500);
}

void DrawAxis() {
  const Vec3f c = {0, 0, 0};
  const float L = 1;
  const Vec3f X = {L, 0, 0}, Y = {0, L, 0}, Z = {0, 0, L};

  glBegin(GL_LINES);
  glColor3f(1, 0, 0);
  glVertex3fv(c.x);
  glVertex3fv((c+X).x);
  glColor3f(0, 1, 0);
  glVertex3fv(c.x);
  glVertex3fv((c+Y).x);
  glColor3f(0, 0, 1);
  glVertex3fv(c.x);
  glVertex3fv((c+Z).x);
  glEnd();
}

bool leftDown = false;
bool rightDown = true;
int origX, origY;

void MouseButton(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_DOWN) {
      leftDown = true;
      rightDown = false;
    } else if (state == GLUT_UP) {
      leftDown = false;
    }
  } else  if (button == GLUT_RIGHT_BUTTON) {
    if (state == GLUT_DOWN) {
      rightDown = true;
    } else if (state == GLUT_UP) {
      rightDown = false;
    }
  }
  if (state == GLUT_DOWN) {
    if (button == 3) {
      eyeDist = eyeDist - 0.04f*eyeDist;
      if (eyeDist <= 0.01f)  eyeDist = 0.01f;
    } else if (button == 4) {
      eyeDist = eyeDist + 0.04f*eyeDist;
    }
  }
  rotateEye(0, 0);
  origX = x;
  origY = y;
  glutPostRedisplay();
}

void MouseMotion(int x, int y) {
  if (leftDown) {
    rotateEye(origX-x, y-origY);
  } else if (rightDown) {
    eyeDist = eyeDist - (y-origY)*0.02f*eyeDist;
    if (eyeDist <= 0.01f)  eyeDist = 0.01f;
    rotateEye(origX-x, 0);
  } else {
    rotateEye(0, 0);
  }
  origX = x;
  origY = y;
  glutPostRedisplay();
}

void Keyboard(unsigned char key, int x, int y) {
  Vec3f forward = eyeDir.unit()*-1;
  Vec3f up = upNorm.unit();
  Vec3f right = (forward^up).unit();
  switch (key) {
    case 'q':
    case 27:  // esc
      exit(0);
      break;
    case 'w':
      center += forward;
      break;
    case 'a':
      center -= right;
      break;
    case 's':
      center -= forward;
      break;
    case 'd':
      center += right;
      break;
    case 'r':
      center -= up;
      break;
    case 'f':
      center += up;
      break;
  }
  rotateEye(0, 0);
  glutPostRedisplay();
}

void SpecialKeyboard(int key, int x, int y) {
  switch (key) {
    case GLUT_KEY_UP:
      rotateEye(0, -1);
      break;
    case GLUT_KEY_DOWN:
      rotateEye(0, 1);
      break;
    case GLUT_KEY_LEFT:
      rotateEye(1, 0);
      break;
    case GLUT_KEY_RIGHT:
      rotateEye(-1, 0);
      break;
  }
  glutPostRedisplay();
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << endl;
    cout << "Usage: ./viewer (filename.obj | -s) [-l]" << endl;
    cout << endl;
    cout << "To load data/test.obj: ./viewer data/test.obj" << endl;
    cout << "To load a custom scene: ./viewer -s" << endl;
    cout << endl;
    return 0;
  }

  // Initialize GLUT
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(window_width, window_height);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("Object viewer");
  glutMouseFunc(MouseButton);
  glutMotionFunc(MouseMotion);
  glutKeyboardFunc(Keyboard);
  glutSpecialFunc(SpecialKeyboard);
  glutDisplayFunc(Display);

  Init();

  if (string(argv[1]) == "-s") {
    char scene_file[] = "obj-data/a/scene.obj";
    cout << "Create scene" << endl;
    disp_scene = true;
    scene_lighting = false;
    initLookup();
    rotateEye(0, 1);
    ParseObj(scene_file, mesh);
    mesh.compute_normals();

    texture_ids = new GLuint[mesh.num_materials()];
    glGenTextures(mesh.num_materials(), texture_ids);

    for (int i = 0; i < mesh.num_materials(); ++i) {
      Material& material = mesh.material(i);
      material.LoadTexture(texture_ids[i]);
    }
  } else if (string(argv[1]) == "-o") {
    char scene_file[] = "obj-data/stargate/original.obj";
    cout << "Custom Object" << endl;
    disp_stargate = true;
    scene_lighting = false;

    ParseObj(scene_file, mesh);
    mesh.compute_normals();

    texture_ids = new GLuint[mesh.num_materials()];
    glGenTextures(mesh.num_materials(), texture_ids);

    for (int i = 0; i < mesh.num_materials(); ++i) {
      Material& material = mesh.material(i);
      material.LoadTexture(texture_ids[i]);
    }
  } else {
    string filename(argv[1]);
    cout << filename << endl;

    // Detect whether to fix the light source(s) to the camera or the world
    scene_lighting = false;
    if (argc > 2 && string(argv[2]) == "-l") {
      scene_lighting = true;
    }

    // Parse the obj file, compute the normals, read the textures

    ParseObj(filename, mesh);
    mesh.compute_normals();
    adjustView();

    texture_ids = new GLuint[mesh.num_materials()];
    glGenTextures(mesh.num_materials(), texture_ids);

    for (int i = 0; i < mesh.num_materials(); ++i) {
      Material& material = mesh.material(i);
      material.LoadTexture(texture_ids[i]);
    }
  }

  glutMainLoop();

  return 0;
}
