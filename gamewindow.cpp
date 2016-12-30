#include "gamewindow.h"

#include <QtGui/QGuiApplication>
#include <QtGui/QMatrix4x4>
//#include <QtGui/QOpenGLShaderProgram>
#include <QtGui/QScreen>

#include <QtCore/qmath.h>
#include <QMouseEvent>
#include <QKeyEvent>
#include <time.h>
#include <sys/time.h>
#include <iostream>

#include <QtCore>
#include <QtGui>

#include <omp.h>
#include "Perlin_Noise/PerlinNoise.cpp"

int numParticules = 5000;
int minP = 0;
int maxP = 360;


using namespace std;

static const char *vertexShaderSource =
  "attribute highp vec4 posAttr;\n"
  "attribute vec3 normal;\n"
  "attribute vec2 texAttr;\n"
  "varying vec2 tex_c;\n"
  "varying vec3 col_normals;\n"
  "uniform lowp vec4 colAttr;\n"
  "varying lowp vec4 col;\n"
  "uniform highp mat4 matrix;\n"
  "varying float height;\n"

  "float abs(float a)\n{\n"
  "if(a < 0.0)\n {a *= -1.0;}\n"
  "return a;\n}\n"
  
  "void main() {\n"
  "float alt = abs(posAttr[2]) * 10.0;\n"
  " col = vec4(alt,0.8,0.3,1.0);\n"
  " tex_c = texAttr;\n"
  " col_normals = normal;\n"
  " gl_Position = matrix * posAttr;\n"
  " height = posAttr[2];\n"
  
  "}\n";

static const char *fragmentShaderSource =
  "varying vec3 col_normals;\n"
  "varying vec2 tex_c;\n"
  "varying lowp vec4 col;\n"
  "uniform sampler2D grass;\n"
  "uniform sampler2D mtn;\n"
  "varying float height;\n"
  
  "void main() {\n"
  "vec4 couleur = vec4(1.0,0.9,0.9,1.0);\n"
  "vec3 direction = vec3(-0.9,0.0,-1.0);\n"
  "vec4 tmp = vec4(0.15,0.1,0.3,1.0);\n"
  "float intens = max(0.0, dot(normalize(col_normals), -direction));\n"
  "vec4 tex_col;\n"
  
  "if(height > 0.160)\n"
  "tex_col = texture2D(mtn, tex_c.st);\n"
  "else\n"
  "{\n"
  " if(height < 0.160 && height > 0.128)\n"
  "{\n"
  "float a = 32.0 - (0.160 - height)*1000.0;\n"
  "float b = (32.0 - a)/32.0;\n"
  "a /= 32.0;\n"
  "tex_col = b*texture2D(grass, tex_c.st) + a*texture2D(mtn, tex_c.st);\n"
  "}\n"
  "else\n"
  "tex_col = texture2D(grass, tex_c.st);\n"
  "}\n"
  "gl_FragColor = col*tex_col * vec4(couleur.rgb * (intens + 0.2), 1.0);\n"
  //"gl_FragColor = tex_col;"
  "}\n";



GLuint GameWindow::loadShader(GLenum type, const char *source)
{
  GLuint shader = glCreateShader(type);
  glShaderSource(shader, 1, &source, 0);
  glCompileShader(shader);
  return shader;
}


GameWindow::GameWindow()
{
  nbTick = 0;
  m_frame = 0;
  master = true;
  
  timer = new QTimer();
  timer->connect(timer, SIGNAL(timeout()),this, SLOT(renderNow()));
  timer->start(maj);
}

void GameWindow::initialize()
{
  init_terrain();
  init_marker();
  init_textures();
  
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glEnable(GL_DEPTH_TEST);
}

void GameWindow::init_textures()
{
  texture = new QOpenGLTexture(QImage("text.jpg").mirrored());
  texture->setMinificationFilter(QOpenGLTexture::Linear);
  texture->setMagnificationFilter(QOpenGLTexture::Linear);
  mountain = new QOpenGLTexture(QImage("grock.jpg").mirrored());
  mountain->setMinificationFilter(QOpenGLTexture::Linear);
  mountain->setMagnificationFilter(QOpenGLTexture::Linear);
  
}

void GameWindow::init_terrain_shader() 
{
  m_program = new QOpenGLShaderProgram(this);
  m_program->addShaderFromSourceCode(QOpenGLShader::Vertex, vertexShaderSource);
  m_program->addShaderFromSourceCode(QOpenGLShader::Fragment, fragmentShaderSource);
  m_program->link();
  m_posAttr = m_program->attributeLocation("posAttr");
  m_colAttr = m_program->uniformLocation("colAttr");
  m_normals = m_program->attributeLocation("normal");
  m_tex = m_program->attributeLocation("texAttr");
  m_matrixUniform = m_program->uniformLocation("matrix");
  m_grass = m_program->uniformLocation("grass");
  m_mtn = m_program->uniformLocation("mtn");
} 

void GameWindow::init_marker_shader() 
{

}

void GameWindow::init_terrain()
{
  init_matrix();
  calc_point(":/heightmap-1.png");
  calc_normals();
  calc_tex();
  
  triangle_terrain = new GLfloat[m_image.width()*2*(m_image.height() - 1)*3];
  
  calc_triangle();
  calc_humid();
  init_terrain_shader();
  
}

void GameWindow::init_marker()
{
  marker_x = 0;
  marker_y = 0;
  
  marker = new GLfloat[6*3*3];
  
  init_point_marker();
}

void GameWindow::init_point_marker()
{
  GLfloat tmp_marker[18*3] = {0.0, 0.0, 0.0,
			      -0.01, 0.01, 0.03,
			      0.01, 0.01, 0.03,
			      0.0, 0.0, 0.0,
			      -0.01, 0.01, 0.03,
			      -0.01, -0.01, 0.03,
			      0.0, 0.0, 0.0,
			      -0.01, -0.01, 0.03,
			      0.01, -0.01, 0.03,
			      0.0, 0.0, 0.0,
			      0.01, -0.01, 0.03,
			      0.01, 0.01, 0.03,
			      0.01, 0.01, 0.03,
			      -0.01, -0.01, 0.03,
			      0.01, -0.01, 0.03,
			      0.01, 0.01, 0.03,
			      -0.01, -0.01, 0.03,
			      -0.01, 0.01, 0.03};
  
  mark_norm = new GLfloat[6*3*3];
  for(int i = 0; i < 18*3; i++)
    {
      mark_norm[i] = 1.0;
      marker[i] = tmp_marker[i];
    }
}

void GameWindow::calc_point(QString localPath)
{
  if (QFile::exists(localPath))
    {
      m_image = QImage(localPath);
      
    }
  unsigned int seed = rand()%255, height = m_image.height(), width = m_image.width();
  PerlinNoise pn(seed);
  

  
  uint id = 0;
  p = new point[width * height];
  QRgb pixel;
  for(int i = 0; i < width; i++)
    {
      for(int j = 0; j < height; j++)
        {
	  double x = (double)j/((double)width);
	  double y = (double)i/((double)height);
	  
	  //  double n = pn.noise(10 * x, 10 * y, 0.8);
	  double n;
	  n =  pn.noise(10 * x, 10 * y, 0.8);
	  //n = 20 * pn.noise(x, y, 0.8);
	  //n = n - floor(n);
	  
	  id = i*width +j;

	  p[id].x = (float)i/(width) - ((float)width/2.0)/width;
	  p[id].y = (float)j/height - ((float)height/2.0)/height;
	  p[id].z = floor(255*n);
	  //if(p[id].z < 128)
	  //p[id].z = 128;
	  p[id].z *= 0.001f;
	  //cout << "z = " << p[id].z << endl;
        }
    }
}

void GameWindow::calc_humid()
{
  unsigned int seed = rand()%255;
  PerlinNoise pn(seed);
  GLfloat points[m_image.width()*m_image.height()];
  int ip, id;
  
  for(int i = 0; i < m_image.width(); i++)
    {
      for(int j = 0; j < m_image.height(); j++)
        {
	  double x = (double)j/((double)m_image.width());
	  double y = (double)i/((double)m_image.height());
	   
	  double n;
	  n =  pn.noise(10 * x, 10 * y, 0.8);
	  points[i * m_image.width() + j] = floor(255*n) * 0.001f;
	}
    }

  humid = new GLfloat[m_image.width()*2*(m_image.height() - 1)];
  
  for(int i = 0; i < m_image.height()-1; i++)
    for(int j = 0; j < m_image.width()*2; j++)
      {
	ip = (j%2 == 1?i+1:i);
	id = ip*m_image.width() +j/2;
	 
	humid[i*m_image.width()*2+j] = points[id];
      }
}
  

void GameWindow::render()
{
  const qreal retinaScale = devicePixelRatio();
  glViewport(0, 0, width() * retinaScale, height() * retinaScale);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  
  m_frame++;

  transform_matrix();
  draw_terrain();
  draw_marker();

}

bool GameWindow::event(QEvent *event)
{
  switch (event->type())
    {
    case QEvent::UpdateRequest:
      renderNow();
      return true;
    default:
      return QWindow::event(event);
    }
}

void GameWindow::init_matrix()
{
  matrix_terrain.ortho(-1.0, 1.0, -1.0, 1.0, -100.0, 100.0);
}

void GameWindow::keyPressEvent(QKeyEvent *event)
{
  int rotation;
  switch(event->key())
    {
    case 'T':
      c->zm += 0.1f;
      break;
    case 'R':
      c->zm -= 0.1f;
      break;
    case 'Z':
      c->ss += 0.10f;
      break;
    case 'S':
      c->ss -= 0.10f;
      break;
    case 'A':
      c->rotX += 1.0f;
      rotation = ((int)(c->rotX)+180) % 361; 
      rotation -= 180;
      c->rotX = (float)rotation;
      cout << "rt = " << c->rotX << endl;
      break;
    case 'E':
      c->rotX -= 1.0f;
      break;
    case 'Q':
      c->rotY += 0.10f;
      break;
    case 'D':
      c->rotY -= 0.10f;
      break;
    case 'O':
      if(marker_y < m_image.height() - 1)
	marker_y++;
      break;
    case 'L' :
      if(marker_y > 0)
	marker_y--;
      break;
    case 'M' :
      if(marker_x < m_image.width() - 1)
	marker_x++;
      break;
    case 'K' :
      if(marker_x > 0)
	marker_x--;
      break;
    case 'W' :
      explosionCrater(idMarker, 20.0f, 0.05f, 0.005f, 0.4f, 0.2f);
      break;
    }
 
}

void GameWindow::transform_matrix()
{
  QMatrix4x4 tf;

  tf.rotate(0.0f,1.0f,0.0f,0.0f);
  tf.rotate(c->rotX,0.0f,0.0f,1.0f);
  tf.translate(c->rotY, c->ss, 0.0f);
  tf.scale(c->zm, c->zm, c->zm);

  tf_terrain = tf;
  
  //tf.setToIdentity();

  tf.translate(0.0f, 0.0f, mr_hover);
  tf.translate((float)marker_y/(float)m_image.width() - 0.5, (float)marker_x/(float)m_image.height() - 0.5, 0.4);
  tf.scale(4.0, 4.0, 4.0);

  tf_marker = tf;
}
  
void GameWindow::draw_marker()
{
  mr_rotat += 0.05f;
  
  if(mr_hover < 0.0)
    hover_s = 1.0;
  if(mr_hover > 0.03)
    hover_s = -1.0;
  
  mr_hover += 0.002f * hover_s;
  
  idMarker = marker_y*m_image.width() + marker_x;
  
  m_program->bind();

  QMatrix4x4 mat = matrix_terrain;
  mat *= tf_marker;
  m_program->setUniformValue(m_matrixUniform, mat);
  
  GLfloat marktex[18*2];
  for (int i = 0; i < 18*2;i++)
    marktex[i] = tex_cord[i];

  glEnableVertexAttribArray(m_posAttr);
  glVertexAttribPointer(m_posAttr, 3, GL_FLOAT, GL_FALSE, 0, marker);
  glEnableVertexAttribArray(m_normals);
  glVertexAttribPointer(m_normals, 3, GL_FLOAT, GL_FALSE, 0, mark_norm);
  glEnableVertexAttribArray(m_normals);
  glVertexAttribPointer(m_tex, 2, GL_FLOAT, GL_FALSE, 0, marktex);

  glUniform1i(m_grass, 0);
  glUniform1i(m_mtn, 1);
  glActiveTexture(GL_TEXTURE0);
  texture->bind();
  glActiveTexture(GL_TEXTURE1);
  mountain->bind();

  glDrawArrays(GL_TRIANGLES, 0, 18);

  glDisableVertexAttribArray(m_posAttr);
  glDisableVertexAttribArray(m_normals);
  glDisableVertexAttribArray(m_tex);  
    
  m_program->release();
}

void GameWindow::explosionCrater(int id, float R, float D, float h, float S, float F)
{
  float n0, n1, m0, m1, w, a, b, c, d, delta;

  for (int x = marker_x -R-1; x < marker_x + R+2; x++){
    for (int y = marker_y -R-1; y< marker_y + R+2; y++){

      int idPoint = y * m_image.width() + x;
      if(idPoint < 0)
	cout << "point = " <<  idPoint << endl;

      float r = sqrt(pow(p[idPoint].x-p[id].x,2)+pow(p[idPoint].y-p[id].y,2))*256;
      float Rn = 2*r/R;

      if (Rn <= 1-F){
	n0 = -1;
	m0 = 0;
	n1 = 0;
	m1 = S*(1-F);
	w = Rn/(1-F);
      } else if (Rn <= (1-F)/2){
	n1 = h/D;
	n0 = 0;
	m1 = 0;
	m0 = S*(F/2)*D/h;
	w = (Rn-(1-F))/(F/2);
      } else if (Rn <= 1){
	n0 = h/D;
	n1 = 0;
	m0 = 0;
	m1 = 0;
	w = (Rn - (1-F/2))/(F/2);
      } else {
	n0 = 0;
	n1 = 0;
	m0 = 0;
	m1 = 0;
	w = 0;
      }

      a = m1 + m0 + 2*(n0-n1);
      b = 3*(n1-n0)-m1-2*m0;
      c = m0;
      d = n0;

      delta = (a*pow(w,3)+b*pow(w,2)+c*w+d)*D;

      //cout << delta << endl;

      p[idPoint].z += delta;
    }
  }
  calc_triangle();
}


void GameWindow::calc_normals()
{
  uint p1 ,p2 ,p3 ,id, ip, *tri = new uint[m_image.width()*2*(m_image.height() - 1)];
  float x, y, z, a, b, c;

  //cout << "done" << endl;

  for(int i = 0; i < m_image.height()-1; i++)
    for(int j = 0; j < m_image.width()*2; j++)
      {
	ip = (j%2 == 1?i+1:i);
	id = ip*m_image.width() +j/2;
	 
	tri[(i*m_image.width()*2+j)] = id;
      }
  
  uint *truetri = new uint[(m_image.width()*2 - 2)*3*(m_image.height() - 1)];
  
  for(int i = 0; i < m_image.height()-1; i++)
    for(int j = 0; j < m_image.width()*2; j++)
      {
	if(!j)
	  {
	    truetri[i*(m_image.width()*2 - 2)*3] = tri[i*m_image.width()*2];
	    id = tri[i*m_image.width()*2 + 1];
	    ip = tri[i*m_image.width()*2 + 2];
 
	    truetri[3*(i*(m_image.width()*2 - 2)) + 1] = id;
	    truetri[3*(i*(m_image.width()*2 - 2)) + 2] = ip;
	  }
	else
	  {
	    truetri[3*(i*(m_image.width()*2 - 2) + j)] = id;
	    truetri[3*(i*(m_image.width()*2 - 2) + j) + 1] = ip;
	    
	    id = ip;
	    ip = tri[i*m_image.width()*2 + j];
	    
	    truetri[3*(i*(m_image.width()*2 - 2) + j) + 2] = ip;
	  }
      }
  
  delete [] tri;

  GLfloat *p_norm = new GLfloat[m_image.width()*m_image.height()*3];
  for(int i = 0; i < m_image.width()*m_image.height()*3; i++)
    p_norm[i] = 0;
  
  for(int i = 0; i < m_image.height()-1; i++)
    for(int j = 0; j < m_image.width()*2 - 2; j++)
      {
	p1 = truetri[3*(i*(m_image.width()*2 - 2) + j)];
	p2 = truetri[3*(i*(m_image.width()*2 - 2) + j) + 1];
	p3 = truetri[3*(i*(m_image.width()*2 - 2) + j) + 2];
  
	x = p[p2].x - p[p1].x;
	y = p[p2].x - p[p1].y;
	z = p[p2].x - p[p1].z;

	a = p[p3].x - p[p1].x;
	b = p[p3].x - p[p1].y;
	c = p[p3].x - p[p1].z;
	
	//cout << " x = " << x << " y = " << y << " z = " << z << " a = " << a << " b = " << b << " c = " << c << endl;
	//cout << " x = " << z*b - y*c << " y = " << z*a - x*c << " z = " << y*a - x*b << endl;

	p_norm[3*p1] += z*b - y*c;
	p_norm[3*p1 + 1] += z*a - x*c;
	p_norm[3*p1 + 2] += y*a - x*b;

	p_norm[3*p3] += z*b - y*c;
	p_norm[3*p3 + 1] += z*a - x*c;
	p_norm[3*p3 + 2] += y*a - x*b;

	p_norm[3*p2] += z*b - y*c;
	p_norm[3*p2 + 1] += z*a - x*c;
	p_norm[3*p2 + 2] += y*a - x*b;
      }
  
  delete[] truetri;

  normals = new GLfloat[m_image.width()*2*(m_image.height() - 1)*3];
  
  for(int ind = 0; ind < m_image.width()*2*(m_image.height() - 1)*3;ind++)
    {
      normals[ind] = 1.0;
    }
  
  
  for(int i = 0; i < m_image.height()-1; i++)
    {
      for(int j = 0; j < m_image.width()*2; j++)
        {
	  ip = (j%2 == 1?i+1:i);
	  id = ip*m_image.width() +j/2;
	  normals[3*(i*m_image.width()*2+j)] = p_norm[3*id];
	  normals[3*(i*m_image.width()*2+j) + 1] = p_norm[3*id + 1];
	  normals[3*(i*m_image.width()*2+j) + 2] = p_norm[3*id + 2];
	}
    }
  
  delete [] p_norm;
  
}

void GameWindow::calc_tex()
{
  
  tex_cord = new GLfloat[m_image.width()*2*(m_image.height() - 1)*2];

  float u = 0.0 ,v = 0.0;
  for(int i = 0; i < m_image.height()-1; i++)
    {
      for(int j = 0; j < m_image.width()*2; j++)
        {
	  tex_cord[2*(i*m_image.width()*2+j)] = u; 
	  tex_cord[2*(i*m_image.width()*2+j) + 1] = v; 
	  
	  if(v == 0.0)
	    {
	      if(u == 0.0)
		v = 1.0;
	      else
		{
		  v = 0.0;
		  u = 0.0;
		}
	    }
	  else
	    {
	      if(u == 0.0)
		{
		  u = 1.0;
		}
	      else
		v = 0.0;		  
	    }
	}
    }
  
}

void GameWindow::calc_triangle()
{
  uint ip,id = 0;
		   
  for(int i = 0; i < m_image.height()-1; i++)
    {
      for(int j = 0; j < m_image.width()*2; j++)
        {
	  
	  ip = (j%2 == 1?i+1:i);
	  id = ip*m_image.width() +j/2;
	  //cout << "i = " <<  p[id].x << " ip = " << p[id].y << " j = " <<"  --- "  << "  width " << p[id+m_image.width()].x<< "  height " << p[id+m_image.height()].y << endl;
	  triangle_terrain[3*(i*m_image.width()*2+j)] = p[id].x;
	  triangle_terrain[3*(i*m_image.width()*2+j) + 1] = p[id].y;
	  triangle_terrain[3*(i*m_image.width()*2+j) + 2] = p[id].z;
	  //cout << "norm =   " << normals[3*(i*m_image.width()*2+j)] << endl;
	}
    }
}	

void GameWindow::draw_terrain()
{
  m_program->bind();

  QMatrix4x4 mat = matrix_terrain;
  mat *= tf_terrain;
  
  m_program->setUniformValue(m_matrixUniform, mat);
 
  glEnableVertexAttribArray(m_posAttr);
  glVertexAttribPointer(m_posAttr, 3, GL_FLOAT, GL_FALSE, 0, triangle_terrain);
  glEnableVertexAttribArray(m_normals);
  glVertexAttribPointer(m_normals, 3, GL_FLOAT, GL_FALSE, 0, normals);
  glEnableVertexAttribArray(m_tex);
  glVertexAttribPointer(m_tex, 2, GL_FLOAT, GL_FALSE, 0, tex_cord);
 
  int cpt = 0;
  glUniform1i(m_grass, 0);
  glUniform1i(m_mtn, 1);
  glActiveTexture(GL_TEXTURE0);
  texture->bind();
  glActiveTexture(GL_TEXTURE1);
  mountain->bind();
  for(cpt = 0;cpt < m_image.height() - 1; cpt++)
    {
      glDrawArrays(GL_TRIANGLE_STRIP, cpt * m_image.width() * 2, m_image.width() * 2);
    }
  
  glDisableVertexAttribArray(m_posAttr);
  glDisableVertexAttribArray(m_normals);
  glDisableVertexAttribArray(m_tex);
    
  m_program->release();
}

