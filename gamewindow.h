#ifndef GAMEWINDOW_H
#define GAMEWINDOW_H
#include <QtGui/QOpenGLShaderProgram>
#include "openglwindow.h"
#include <QOpenGLTexture>

struct point
{
    float x, y ,z;
};

class paramCamera
{
public:
    float rotX = -45.0;
    float rotY = -45.0;
    float ss = 1.0f;
    float anim = 0.0f;

    int etat = 2;
};

class GameWindow : public OpenGLWindow
{
Q_OBJECT

public:
    GameWindow();
    GameWindow(int maj);
    void initialize();
    void render();
    bool event(QEvent *event);

    void keyPressEvent(QKeyEvent *event);

    void displayTriangles();
    void displayLines();
    void displayTrianglesC();
    void displayPoints();
    void displayTrianglesTexture();

    void displayColor(float);

    void loadMap(QString localPath);
    void updateParticlesAut();
    void updateParticlesHiv();
    paramCamera* c;

    void setSeason(int );

    void calc_normals();
    void calc_tex();

    GLuint loadShader(GLenum type, const char *source);
public slots:
    void updateSeason();

private:
    int nbTick = 0;
    int m_frame = 0;
    int season, day;
    point* particules;
    bool master = false;

    QImage m_image;
    point *p;
    GLfloat* normals;
    GLfloat* tex;
    int carte=1;
    int maj = 20;

    QOpenGLTexture *texture;

    QTimer *timer;
    QTimer *timerFPS;
    
    GLuint m_posAttr;
    GLuint m_colAttr;
    GLuint m_matrixUniform;
    GLuint m_normals;
    GLuint m_tex;
    
    QOpenGLShaderProgram *m_program;
    
};



#endif // GameWindow_H
