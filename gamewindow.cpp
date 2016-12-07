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
  "attribute vec4 normals;\n"
  "varying vec3 col_normals;\n"
  "uniform lowp vec4 colAttr;\n"
  "varying lowp vec4 col;\n"
  "uniform highp mat4 matrix;\n"
  "void main() {\n"
  /*
    "tmp[1] = posAttr[2]*25.5;"
    "if(tmp[1] < 0.0)"
    "  tmp[1] = 0.0;"
    "if(tmp[1] > 1.0)"
    "  tmp[1] = 1.0;"
    " col = tmp;\n"
  */
  
  " col_normals = normals.xyz;\n"
  " gl_Position = matrix * posAttr;\n"
  "}\n";

static const char *fragmentShaderSource =
  "varying vec3 col_normals;\n"
  "varying lowp vec4 col;\n"
  "void main() {\n"
  "vec3 couleur = vec3(1.0,1.0,1.0);\n"
  "vec3 direction = vec3(-0.5,-0.5,-0.5);\n"
  "vec4 tmp = vec4(0.15,0.1,0.3,1);\n"
  "float intens = max(0.0, dot(normalize(col_normals), -direction));\n"
  "gl_FragColor = tmp*vec4(couleur * (intens + 0.2), 1.0);\n"
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
  maj = 20;
  QString s ("FPS : ");
  s += QString::number(1000/maj);
  s += "(";
  s += QString::number(maj);
  s += ")";
  setTitle(s);
  timer = new QTimer();
  timer->connect(timer, SIGNAL(timeout()),this, SLOT(renderNow()));
  timer->start(maj);

  master = true;
}
GameWindow::GameWindow(int _maj)
{
  nbTick = 0;
  m_frame = 0;
  maj = _maj;
  QString s ("FPS : ");
  s += QString::number(1000/maj);
  s += "(";
  s += QString::number(maj);
  s += ")";
  setTitle(s);
  timer = new QTimer();
  timer->connect(timer, SIGNAL(timeout()),this, SLOT(renderNow()));
  timer->start(maj);
}

void GameWindow::setSeason(int i)
{
  season = i;

  if (i==0) day=80;
  else if (i==1) day = 170;
  else if (i==2) day = 260;
  else if (i==3) day = 350;
}

void GameWindow::updateSeason()
{
  day = (day + 1) % 365;

  if (day==80) season = 0;
  else if (day==170) season = 1;
  else if (day==260) season = 2;
  else if (day==350) season = 3;
}


void GameWindow::initialize()
{
  const qreal retinaScale = devicePixelRatio();

  calc_normals();

  glViewport(0, 0, width() * retinaScale, height() * retinaScale);

  glClearColor(0.0, 0.0, 0.0, 0.0);
  //glMatrixMode(GL_PROJECTION);
  //glLoadIdentity();
  //glOrtho(-1.0, 1.0, -1.0, 1.0, -100.0, 100.0);

  m_program = new QOpenGLShaderProgram(this);
  m_program->addShaderFromSourceCode(QOpenGLShader::Vertex, vertexShaderSource);
  m_program->addShaderFromSourceCode(QOpenGLShader::Fragment, fragmentShaderSource);
  m_program->link();
  m_posAttr = m_program->attributeLocation("posAttr");
  m_colAttr = m_program->uniformLocation("colAttr");
  m_normals = m_program->attributeLocation("normals");
  m_matrixUniform = m_program->uniformLocation("matrix");
  
  
  
  loadMap(":/heightmap-1.png");

  particules = new point[numParticules];

  for(int i = 0; i < numParticules; i++)
    {
      int angle = minP + (rand() % (int)(maxP - minP + 1));
      int dist = (rand() % (int)(100));
      int alt = (rand() % (int)(100));
      float x = dist*sin(
			 ((3.14159 * 2) *
			  angle
			  )/360
			 );
      float y = dist*cos(
			 ((3.14159 * 2) *
			  angle
			  )/360
			 );

      // x and y are in (-100,100)

      particules[i].x = (float)(x)/(m_image.width());
      particules[i].y = (float)(y)/(m_image.height());
      particules[i].z = (float)(alt)/100;
    }

}

///*
void GameWindow::loadMap(QString localPath)
{
  
  if (QFile::exists(localPath)) {
    m_image = QImage(localPath);
  }
  unsigned int seed = 237, height = 257, width = 257;
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
	  if(p[id].z < 128)
	    p[id].z = 128;
	  p[id].z *= 0.001f;
	  //cout << "z = " << p[id].z << endl;
        }
    }
}
//*/
/*
void GameWindow::loadMap(QString localPath)
{

  if (QFile::exists(localPath)) {
    m_image = QImage(localPath);
  }


  uint id = 0;
  p = new point[m_image.width() * m_image.height()];
  QRgb pixel;
  for(int i = 0; i < m_image.width(); i++)
    {
      for(int j = 0; j < m_image.height(); j++)
        {

	  pixel = m_image.pixel(i,j);

	  id = i*m_image.width() +j;

	  p[id].x = (float)i/(m_image.width()) - ((float)m_image.width()/2.0)/m_image.width();
	  p[id].y = (float)j/(m_image.height()) - ((float)m_image.height()/2.0)/m_image.height();
	  p[id].z = 0.001f * (float)(qRed(pixel));

        }
    }
}
*/

void GameWindow::render()
{
  nbTick += maj;

  if(nbTick >= 1000)
    {
      QString s ("FPS : ");
      s += QString::number(m_frame);
      s += "(th=";
      s += QString::number(1000/maj);
      s += "-";
      s += QString::number(maj);
      s += ")";
      s += " day ";
      s += QString::number(day);
      if (season==0) s += " Printemps ";
      else if (season==1) s += " Eté ";
      if (season==2) s += " Automne ";
      if (season==3) s += " Hiver ";
      setTitle(s);
      nbTick = 0;
      m_frame = 0;
    }
  //  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  glClear(GL_COLOR_BUFFER_BIT);
  
  
   glLoadIdentity();
  //glScalef(c->ss,c->ss,c->ss);
  
  glRotatef(c->rotX,1.0f,0.0f,0.0f);
  



  if(c->anim == 0.0f)
    {
      glRotatef(c->rotY,0.0f,0.0f,1.0f);
    }
  else
    {
        glRotatef(c->anim,0.0f,0.0f,1.0f);
      if(master)
	c->anim +=0.05f;
    }



  switch(c->etat)
    {
    case 0:
      displayPoints();
      break;
    case 1:
      displayLines();
      break;
    case 2:
      displayTriangles();
      break;
    case 3:
      displayTrianglesC();
      break;
    case 4:
      displayTrianglesTexture();
      break;
    case 5:

      displayTrianglesTexture();
      displayLines();
      break;
    default:
      displayPoints();
      break;
    }
  if (season == 2)
    updateParticlesAut();
  else if (season == 3)
    updateParticlesHiv();

  m_frame++;

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

void GameWindow::keyPressEvent(QKeyEvent *event)
{
  switch(event->key())
    {

    case 'C':
      if(c->anim == 0.0f)
	c->anim = c->rotY;
      else
	c->anim = 0.0f;
      break;
    case 'Z':
      c->ss += 0.10f;
      break;
    case 'S':
      c->ss -= 0.10f;
      break;
    case 'A':
      c->rotX += 1.0f;
      break;
    case 'E':
      c->rotX -= 1.0f;
      break;
    case 'Q':
      c->rotY += 1.0f;
      break;
    case 'D':
      c->rotY -= 1.0f;
      break;
    case 'W':
      c->etat ++;
      if(c->etat > 5)
	c->etat = 0;
      break;
    case 'P':
      maj++;
      timer->stop();
      timer->start(maj);
      break;
    case 'O':
      maj--;
      if(maj < 1)
	maj = 1;
      timer->stop();
      timer->start(maj);
      break;
    case 'L':
      maj = maj - 20;
      if(maj < 1)
	maj = 1;
      timer->stop();
      timer->start(maj);
      break;
    case 'M':
      maj = maj + 20;

      timer->stop();
      timer->start(maj);
      break;
    case 'X':
      carte ++;
      if(carte > 3)
	carte = 1;
      QString depth (":/heightmap-");
      depth += QString::number(carte) ;
      depth += ".png" ;

      loadMap(depth);
      break;


    }

}


void GameWindow::displayPoints()
{
  m_program->bind();
  QMatrix4x4 matrix;
  //  matrix.ortho(-1.0, 1.0, -1.0, 1.0, -100.0, 100.0);
  matrix.rotate(c->rotX,1.0f,0.0f,0.0f);
  matrix.rotate(c->anim,0.0f,0.0f,1.0f);
  
  m_program->setUniformValue(m_matrixUniform, matrix);
  QColor cl(20,10,10,55);
  m_program->setUniformValue(m_colAttr,cl);

  

  GLfloat vertices[m_image.width()*m_image.height()*3];
  
  uint id = 0;
  for(int i = 0; i < m_image.width(); i++)
    {
      for(int j = 0; j < m_image.height(); j++)
	{
	 id = i*m_image.width() +j;
	 vertices[3*id] = p[id].x;
	 vertices[3*id + 1] = p[id].y;
	 vertices[3*id + 2] = p[id].z;
	}
    }
  
  GLfloat colors[] = {
    1.0f, 0.0f, 1.0f,
      };
       
  glVertexAttribPointer(m_posAttr, 3, GL_FLOAT, GL_FALSE, 0, vertices);
  //  glVertexAttribPointer(m_colAttr, 3, GL_FLOAT, GL_FALSE, 0, colors);
  glVertexAttribPointer(m_normals, 3, GL_FLOAT, GL_FALSE, 0, normals);
     
  glEnableVertexAttribArray(0);
  //EnableVertexAttribArray(1);
  
  
  
  glDrawArrays(GL_POINTS, 0, m_image.width()*m_image.height());
 
  //DisableVertexAttribArray(1);
  glDisableVertexAttribArray(0);
    
  m_program->release();
 
}

void GameWindow::calc_normals()
{
  normals = new GLfloat[m_image.width()*m_image.height()*3];
  for(int i = 0; i < m_image.width()*m_image.height()*3;i++)
    normals[i] = 0;

  uint p1 ,p2 ,p3 ,id, ip, tri[m_image.width()*2*(m_image.height() - 1)];
  float x, y, z, a, b, c;

  //cout << "done" << endl;

  for(int i = 0; i < m_image.height()-1; i++)
    for(int j = 0; j < m_image.width()*2; j++)
      {
	ip = (j%2 == 1?i+1:i);
	id = ip*m_image.width() +j/2;
	 
	tri[(i*m_image.width()*2+j)] = id;
      }
  
  uint truetri[(m_image.width()*2 - 2)*3*(m_image.height() - 1)];
  
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

  for(int i = 0; i < m_image.height()-1; i++)
    for(int j = 0; j < m_image.width()*2 - 2; j++)
      {
	p1 = truetri[3*(i*(m_image.width()*2 - 2) + j)];
	p2 = truetri[3*(i*(m_image.width()*2 - 2) + j) + 1];
	p3 = truetri[3*(i*(m_image.width()*2 - 2) + j) + 2];
  
	x = p[p2].x - p[p1].x;
	y = p[p2].x - p[p1].y;
	z = p[p2].x - p[p1].z;

	a = p[p2].x - p[p1].x;
	b = p[p2].x - p[p1].y;
	c = p[p2].x - p[p1].z;
	
	normals[3*p1] += z*b - y*c;
	normals[3*p1 + 1] += z*a - x*c;
	normals[3*p1 + 2] += y*a - x*b;
      }
}
	
void GameWindow::displayTriangles()
{
 
  
   
  //m_program->setUniformValue(m_posAttr,cl);
  
  

  GLfloat vertices[m_image.width()*2*(m_image.height() - 1)*3];
  uint ip,id = 0;
		   
  for(int i = 0; i < m_image.height()-1; i++)
    {
      for(int j = 0; j < m_image.width()*2; j++)
        {
	  
	  ip = (j%2 == 1?i+1:i);
	  id = ip*m_image.width() +j/2;
	  //cout << "i = " <<  p[id].x << " ip = " << p[id].y << " j = " <<"  --- "  << "  width " << p[id+m_image.width()].x<< "  height " << p[id+m_image.height()].y << endl;
	  vertices[3*(i*m_image.width()*2+j)] = p[id].x;
	  vertices[3*(i*m_image.width()*2+j) + 1] = p[id].y;
	  vertices[3*(i*m_image.width()*2+j) + 2] = p[id].z;
	}
    }

  m_program->bind();
  QMatrix4x4 matrix;
  //  matrix.ortho(-1.0, 1.0, -1.0, 1.0, -100.0, 100.0);
  matrix.rotate(c->rotX,1.0f,0.0f,0.0f);
  matrix.rotate(c->anim,0.0f,0.0f,1.0f);
  
  m_program->setUniformValue(m_matrixUniform, matrix);
  
  
  GLfloat colors[] = {
        1.0f, 1.0f, 0.0f,
      };
       
  glVertexAttribPointer(m_posAttr, 3, GL_FLOAT, GL_FALSE, 0, vertices);
  glVertexAttribPointer(m_normals, 3, GL_FLOAT, GL_FALSE, 0, normals);
       
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);

  /*
  QOpenGLTexture *texture = new QOpenGLTexture(QImage("text.jpg").mirrored());
  texture->setMinificationFilter(QOpenGLTexture::LinearMipMapLinear);
  texture->setMagnificationFilter(QOpenGLTexture::Linear);
  */

  int cpt = 0;
  //  texture->bind();
  for(;cpt < m_image.height() - 1; cpt++)
    {
      glDrawArrays(GL_TRIANGLE_STRIP, cpt * m_image.width() * 2, m_image.width() * 2);
    }
  
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(0);
    
  m_program->release();
  


  /*
  id = i*m_image.width() +(j+1);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = (i+1)*m_image.width() +j;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);



	  id = i*m_image.width() +(j+1);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = (i+1)*m_image.width() +j+1;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = (i+1)*m_image.width() +j;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
        }
    }

  glEnd();
  */
}

void GameWindow::displayTrianglesC()
{
  glColor3f(1.0f, 1.0f, 1.0f);
  glBegin(GL_TRIANGLES);
  uint id = 0;

  for(int i = 0; i < m_image.width()-1; i++)
    {
      for(int j = 0; j < m_image.height()-1; j++)
        {
	  glColor3f(0.0f, 1.0f, 0.0f);
	  id = i*m_image.width() +j;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = i*m_image.width() +(j+1);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = (i+1)*m_image.width() +j;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);


	  glColor3f(1.0f, 1.0f, 1.0f);
	  id = i*m_image.width() +(j+1);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = (i+1)*m_image.width() +j+1;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = (i+1)*m_image.width() +j;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
        }
    }
  glEnd();
}


void GameWindow::displayLines()
{
  //glColor3f(1.0f, 1.0f, 1.0f);
  glLineWidth(0.05);
  // glBegin(GL_LINES);
  m_program->bind();
  QMatrix4x4 matrix;
  //  matrix.ortho(-1.0, 1.0, -1.0, 1.0, -100.0, 100.0);
  matrix.rotate(c->rotX,1.0f,0.0f,0.0f);
  matrix.rotate(c->anim,0.0f,0.0f,1.0f);
  
  m_program->setUniformValue(m_matrixUniform, matrix);
  

  GLfloat vertices[m_image.width()*m_image.height()*3*2];
  
  uint id2=0,id = 0;
  for(int i = 0; i < m_image.width(); i++)
    {
      for(int j = 0; j < m_image.height(); j++)
	{
	 id = i*m_image.width() +j;
	 id2 = j*m_image.height() + i ;
	 vertices[3*id] = p[id].x;
	 vertices[3*id + 1] = p[id].y;
	 vertices[3*id + 2] = p[id].z;
	 
	 vertices[3*id2  + m_image.width()*m_image.height()*3] = p[id2].x;
	 vertices[3*id2 + 1 + m_image.width()*m_image.height()*3] = p[id2].y;
	 vertices[3*id2 + 2 + m_image.width()*m_image.height()*3] = p[id2].z;
	}
    }

    GLfloat colors[] = {
      1.0f, 0.3f, 0.0f
      };
       
  glVertexAttribPointer(m_posAttr, 3, GL_FLOAT, GL_FALSE, 0, vertices);
  glVertexAttribPointer(m_colAttr, 3, GL_FLOAT, GL_FALSE, 0, colors);
  glVertexAttribPointer(m_normals, 3, GL_FLOAT, GL_FALSE, 0, normals);
       
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  
  int cpt = 0;
  for(cpt; cpt < m_image.height();cpt++)
    glDrawArrays(GL_LINE_STRIP, cpt*m_image.width(), m_image.width());
   for(cpt; cpt < m_image.width();cpt++)
     glDrawArrays(GL_LINE_STRIP, cpt*m_image.height(), m_image.height());
  

  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(0);
    
  m_program->release();


  
  /*
  for(int i = 0; i < m_image.width()-1; i++)
    {
      for(int j = 0; j < m_image.height()-1; j++)
        {

	  id = i*m_image.width() +j;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = i*m_image.width() +(j+1);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);

	  id = (i+1)*m_image.width() +j;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = i*m_image.width() +j;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);

	  id = (i+1)*m_image.width() +j;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = i*m_image.width() +(j+1);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);

	  id = i*m_image.width() +(j+1);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = (i+1)*m_image.width() +j+1;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);

	  id = (i+1)*m_image.width() +j+1;
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);

	  id = (i+1)*m_image.width() +(j);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
        }
    }

    glEnd();*/
}

void GameWindow::displayTrianglesTexture()
{
  glColor3f(1.0f, 1.0f, 1.0f);
  glBegin(GL_TRIANGLES);
  uint id = 0;

  for(int i = 0; i < m_image.width()-1; i++)
    {
      for(int j = 0; j < m_image.height()-1; j++)
        {

	  id = i*m_image.width() +j;
	  displayColor(p[id].z);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = i*m_image.width() +(j+1);
	  displayColor(p[id].z);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = (i+1)*m_image.width() +j;
	  displayColor(p[id].z);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);



	  id = i*m_image.width() +(j+1);
	  displayColor(p[id].z);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = (i+1)*m_image.width() +j+1;
	  displayColor(p[id].z);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
	  id = (i+1)*m_image.width() +j;
	  displayColor(p[id].z);
	  glVertex3f(
		     p[id].x,
		     p[id].y,
		     p[id].z);
        }
    }
  glEnd();
}


void GameWindow::displayColor(float alt)
{
  if (alt > 0.2)
    {
      glColor3f(01.0f, 1.0f, 1.0f);
    }
  else     if (alt > 0.1)
    {
      glColor3f(alt, 1.0f, 1.0f);
    }
  else     if (alt > 0.05f)
    {
      glColor3f(01.0f, alt, alt);
    }
  else
    {
      glColor3f(0.0f, 0.0f, 1.0f);
    }

}


void GameWindow::updateParticlesAut()
{
  int id2;
#pragma omp parallel
  {
#pragma omp for
    for(int id = 0; id < numParticules; id++)
      {
	particules[id].z -= 0.0003f * ((float) minP + (rand() % (int)(maxP - minP + 1)));
	id2 = m_image.width()*m_image.width()/4 + (particules[id].x)*m_image.width() + particules[id].y;

	if (id2<0)
	  qDebug() << "error x = " << particules[id].x << "  / " << m_image.width() << " y = " << particules[id].y << " / " << m_image.height();

	if (id2>m_image.width()*m_image.height())
	  qDebug() << "error x = " << particules[id].x << "  / " << m_image.width() << " y = " << particules[id].y << " / " << m_image.height();

	// restart when touching the ground
	if(particules[id].z < p[id2].z)
	  {
	    int angle =minP + (rand() % (int)(maxP - minP + 1));
	    int dist = (rand() % (int)(100 ));
	    int alt = (rand() % (int)(100));
	    float x = dist*sin(
			       ((3.14159 * 2) *
				angle
				)/360
			       );
	    float y = dist*cos(
			       ((3.14159 * 2) *
				angle
				)/360
			       );

	    particules[id].x = (float)(x)/(m_image.width());
	    particules[id].y = (float)(y)/(m_image.height());
	    particules[id].z = (float)(alt)/100;

	  }
	// else display the river or cover the round with snow

      }
  }
  m_program->bind();
  QMatrix4x4 matrix;
  //  matrix.ortho(-1.0, 1.0, -1.0, 1.0, -100.0, 100.0);
  matrix.rotate(c->rotX,1.0f,0.0f,0.0f);
  matrix.rotate(c->anim,0.0f,0.0f,1.0f);
  
  m_program->setUniformValue(m_matrixUniform, matrix);
  
  
  GLfloat vertices[numParticules*3];
  // glColor3f(0.2f, 0.2f, 1.0f);
  //glPointSize(0.01f);
  for(int id = 0; id < numParticules; id++)
    {
      vertices[3*id] =  particules[id].x;
      vertices[3*id + 1] = particules[id].y;
      vertices[3*id + 2] = particules[id].z;
    }
    GLfloat colors[] = {
    0.0f, 0.3f, 1.0f,
      };
       
  glVertexAttribPointer(m_posAttr, 3, GL_FLOAT, GL_FALSE, 0, vertices);
  glVertexAttribPointer(m_colAttr, 3, GL_FLOAT, GL_FALSE, 0, colors);
       
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  
  
  
  glDrawArrays(GL_POINTS, 0, numParticules);
 
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(0);
    
  m_program->release();
  
}

void GameWindow::updateParticlesHiv()
{
  int id2;
#pragma omp parallel
  {
#pragma omp for
    for(int id = 0; id < numParticules; id++)
      {
	particules[id].z -= 0.00001f * ((float) minP + (rand() % (int)(maxP - minP + 1)));
	id2 = m_image.width()*m_image.width()/4 + (particules[id].x)*m_image.width() + particules[id].y;

	if (id2<0)
	  qDebug() << "error x = " << particules[id].x << "  / " << m_image.width() << " y = " << particules[id].y << " / " << m_image.height();

	if (id2>m_image.width()*m_image.height())
	  qDebug() << "error x = " << particules[id].x << "  / " << m_image.width() << " y = " << particules[id].y << " / " << m_image.height();

	// restart when touching the ground
	if(particules[id].z < p[id2].z)
	  {
	    int angle =minP + (rand() % (int)(maxP - minP + 1));
	    int dist = (rand() % (int)(100));
	    int alt =  (rand() % (int)(100));
	    float x = dist*sin(
			       ((3.14159 * 2) *
				angle
				)/360
			       );
	    float y = dist*cos(
			       ((3.14159 * 2) *
				angle
				)/360
			       );

	    particules[id].x = (float)(x)/(m_image.width());
	    particules[id].y = (float)(y)/(m_image.height());
	    particules[id].z = (float)(alt)/100;
	  }
      }
  }
  // glColor3f(1.0f, 1.0f, 1.0f);
  //  glPointSize(0.0001f);
   m_program->bind();
  QMatrix4x4 matrix;
  //  matrix.ortho(-1.0, 1.0, -1.0, 1.0, -100.0, 100.0);
  matrix.rotate(c->rotX,1.0f,0.0f,0.0f);
  matrix.rotate(c->anim,0.0f,0.0f,1.0f);
  
  m_program->setUniformValue(m_matrixUniform, matrix);
  
  
  GLfloat vertices[numParticules*3];
  // glColor3f(0.2f, 0.2f, 1.0f);
  //glPointSize(0.01f);
  for(int id = 0; id < numParticules; id++)
    {
      vertices[3*id] =  particules[id].x;
      vertices[3*id + 1] = particules[id].y;
      vertices[3*id + 2] = particules[id].z;
    }
    GLfloat colors[] = {
      1.0f, 1.0f, 1.0f
  };
       
  glVertexAttribPointer(m_posAttr, 3, GL_FLOAT, GL_FALSE, 0, vertices);
  glVertexAttribPointer(m_colAttr, 3, GL_FLOAT, GL_FALSE, 0, colors);
       
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  
  
  
  glDrawArrays(GL_POINTS, 0, numParticules);
 
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(0);
    
  m_program->release();
}
