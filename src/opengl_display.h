#ifndef OPENGL_DISPLAY_H
#define OPENGL_DISPLAY_H

#include <SDL2/SDL.h>
#include <GL/gl.h>
#include <cstdint>

/**
 * OpenGL窗口显示类
 * 用于在OpenGL窗口中实时显示渲染结果
 */
class OpenGLDisplay {
public:
    OpenGLDisplay(int width, int height, const char* title = "BHRaytracer");
    ~OpenGLDisplay();

    // 初始化OpenGL窗口
    bool init();
    
    // 更新显示（将像素数据上传到OpenGL纹理并显示）
    void update(const char* pixels, int width, int height);
    
    // 处理窗口事件（返回false表示窗口关闭）
    bool handleEvents();
    
    // 交换缓冲区
    void swapBuffers();
    
    // 检查窗口是否应该关闭
    bool shouldClose() const { return m_shouldClose; }
    
    // 清理资源
    void cleanup();

private:
    SDL_Window* m_window;
    SDL_GLContext m_context;
    GLuint m_texture;
    int m_width;
    int m_height;
    bool m_shouldClose;
    bool m_initialized;
    
    // 创建OpenGL纹理
    void createTexture();
};

#endif // OPENGL_DISPLAY_H

