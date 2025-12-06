#include "opengl_display.h"
#include <iostream>
#include <cstring>

OpenGLDisplay::OpenGLDisplay(int width, int height, const char* title)
    : m_window(nullptr), m_context(nullptr), m_texture(0),
      m_width(width), m_height(height), m_shouldClose(false), m_initialized(false) {
}

OpenGLDisplay::~OpenGLDisplay() {
    cleanup();
}

bool OpenGLDisplay::init() {
    if (m_initialized) {
        return true;
    }

    // 初始化SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "SDL初始化失败: " << SDL_GetError() << std::endl;
        return false;
    }

    // 设置OpenGL属性
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

    // 创建窗口
    m_window = SDL_CreateWindow(
        "BHRaytracer",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        m_width,
        m_height,
        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE
    );

    if (!m_window) {
        std::cerr << "窗口创建失败: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return false;
    }

    // 创建OpenGL上下文
    m_context = SDL_GL_CreateContext(m_window);
    if (!m_context) {
        std::cerr << "OpenGL上下文创建失败: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(m_window);
        SDL_Quit();
        return false;
    }

    // 初始化OpenGL
    glEnable(GL_TEXTURE_2D);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    
    // 设置视口
    glViewport(0, 0, m_width, m_height);
    
    // 设置2D投影矩阵
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, m_width, m_height, 0, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // 创建纹理
    createTexture();

    m_initialized = true;
    std::cout << "OpenGL窗口初始化成功: " << m_width << "x" << m_height << std::endl;
    return true;
}

void OpenGLDisplay::createTexture() {
    glGenTextures(1, &m_texture);
    glBindTexture(GL_TEXTURE_2D, m_texture);
    
    // 设置纹理参数
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    
    // 分配纹理内存（初始为空）
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, m_width, m_height, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
}

void OpenGLDisplay::update(const char* pixels, int width, int height) {
    if (!m_initialized || !pixels) {
        return;
    }

    // 更新纹理数据
    glBindTexture(GL_TEXTURE_2D, m_texture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

    // 清除屏幕
    glClear(GL_COLOR_BUFFER_BIT);

    // 绘制纹理到全屏
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, m_texture);
    
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f); glVertex2f(0.0f, 0.0f);
    glTexCoord2f(1.0f, 0.0f); glVertex2f(width, 0.0f);
    glTexCoord2f(1.0f, 1.0f); glVertex2f(width, height);
    glTexCoord2f(0.0f, 1.0f); glVertex2f(0.0f, height);
    glEnd();
    
    glDisable(GL_TEXTURE_2D);
}

bool OpenGLDisplay::handleEvents() {
    if (!m_initialized) {
        return true;
    }

    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        switch (event.type) {
            case SDL_QUIT:
                m_shouldClose = true;
                return false;
            case SDL_KEYDOWN:
                if (event.key.keysym.sym == SDLK_ESCAPE || event.key.keysym.sym == SDLK_q) {
                    m_shouldClose = true;
                    return false;
                }
                break;
            case SDL_WINDOWEVENT:
                if (event.window.event == SDL_WINDOWEVENT_RESIZED) {
                    m_width = event.window.data1;
                    m_height = event.window.data2;
                    glViewport(0, 0, m_width, m_height);
                }
                break;
        }
    }
    return true;
}

void OpenGLDisplay::swapBuffers() {
    if (m_initialized && m_window) {
        SDL_GL_SwapWindow(m_window);
    }
}

void OpenGLDisplay::cleanup() {
    if (m_texture != 0) {
        glDeleteTextures(1, &m_texture);
        m_texture = 0;
    }

    if (m_context) {
        SDL_GL_DeleteContext(m_context);
        m_context = nullptr;
    }

    if (m_window) {
        SDL_DestroyWindow(m_window);
        m_window = nullptr;
    }

    if (m_initialized) {
        SDL_Quit();
        m_initialized = false;
    }
}

