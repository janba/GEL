/*
 *  ShadowBuffer.cpp
 *  02564_Framework
 *
 *  Created by J. Andreas Bærentzen on 01/02/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include "ShadowBuffer.h"

using namespace std;

void ShadowBuffer::relinquish()
{
	glDeleteTextures(1, &dtex);
	glDeleteRenderbuffers(1, &rb);
	glDeleteFramebuffers(1, &fbo);
}

void ShadowBuffer::initialize()
{
	relinquish();
	
	glGenTextures(1, &dtex);
	glBindTexture(GL_TEXTURE_2D, dtex);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LESS);
	
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, dim, dim, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, 0);
	
    glGenFramebuffers(1,&fbo);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
	
	
    glGenRenderbuffers(1,&rb);
    glBindRenderbuffer(GL_RENDERBUFFER, rb);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA, dim, dim);
	
    glFramebufferRenderbuffer(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, rb);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, dtex, 0);
	
	
    if(glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
		cout << "Something wrong with FBO" << endl;
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	
}

ShadowBuffer::ShadowBuffer(int _dim): dim(_dim), fbo(0), rb(0), dtex(0)
{
	initialize();
};

void ShadowBuffer::bind_textures(int dtex_unit)
{
	glActiveTexture(GL_TEXTURE0+dtex_unit);
	glBindTexture(GL_TEXTURE_2D, dtex);
}

GLint ShadowBuffer::enable()
{
	GLint old_draw_buffer;
	glGetIntegerv(GL_DRAW_BUFFER, &old_draw_buffer);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
	
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	return old_draw_buffer;
}
