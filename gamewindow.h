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

    void draw_terrain();
    void draw_marker();
 
    paramCamera* c;

    void calc_point(QString localPath);
    void calc_triangle();
    void calc_normals();
    void calc_tex();
    void calc_humid();
    
    void init_terrain_shader();
    void init_marker_shader();

    void init_terrain();
    void init_marker();

    void init_textures();

    void init_matrix();
    void transform_matrix();

    void init_point_marker();
    void explosionCrater(int, float, float, float, float, float);

    GLuint loadShader(GLenum type, const char *source);

private:
    int nbTick = 0;
    int m_frame = 0;
    
    int maj = 20;

    QTimer *timer;
    QTimer *timerFPS;

    QImage m_image;
    point *p;

    GLfloat* triangle_terrain;
    GLfloat* normals;
    GLfloat* tex_cord;
    GLfloat* humid;

    bool master;

    QOpenGLTexture *texture;
    QOpenGLTexture *mountain;

    int idMarker;
    int marker_x;
    int marker_y;

    GLfloat* marker;
    GLfloat* mark_norm;
    
    QMatrix4x4 matrix_terrain;
    QMatrix4x4 matrix_marker;
    QMatrix4x4 tf_terrain;
    QMatrix4x4 tf_marker;    
    
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
