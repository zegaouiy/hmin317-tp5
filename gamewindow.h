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
    float rotY = 0.0f;
    float ss = 0.0f;
    float anim = 0.0f;
    float zm = 1.0f;

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
    void calc_humid();

    void initMarker();
    void displayExplosionMarker(int);
    void explosionCrater(int, float, float, float, float, float);

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
    GLfloat* tex_cord;
    GLfloat* humid;
    int carte=1;
    int maj = 20;

    QOpenGLTexture *texture;
    QOpenGLTexture *mountain;

    QTimer *timer;
    QTimer *timerFPS;
    
    int idMarker;
    int marker_x;
    int marker_y;
    GLfloat* marker;
    GLfloat* mark_norm;
    GLfloat mr_rotat = 0;
    GLfloat mr_hover = 0;
    GLfloat hover_s = 1.0;

    GLuint m_posAttr;
    GLuint m_colAttr;
    GLuint m_matrixUniform;
    GLuint m_normals;
    GLuint m_tex;
    GLuint m_grass;
    GLuint m_mtn;
    GLuint m_humid;
    
    
    QOpenGLShaderProgram *m_program;
    
};



#endif // GameWindow_H
