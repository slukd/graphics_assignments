#include "MovableGLM.h"
#include <stdio.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>


static void printMat(const glm::mat4 mat)
{
	printf(" matrix: \n");
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			printf("%f ", mat[j][i]);
		printf("\n");
	}
}

MovableGLM::MovableGLM()
{
	ZeroTrans();
}

glm::mat4 MovableGLM::MakeTrans(glm::mat4 &prevTransformations) const
{
	return prevTransformations * MakeTrans();
}

glm::mat4 MovableGLM::MakeTrans() const
{
	return   trans * rot * scl ;
}

void MovableGLM::MyTranslate(glm::vec3 delta,int mode)
{
	trans = glm::translate(trans,delta);
}

void  MovableGLM::MyRotate(float angle,glm::vec3 &vec,int mode)
{
	if (mode == 1) {
		glm::mat4 inverse_rot = glm::inverse(get_rot());
		glm::vec3 new_vec = glm::normalize(glm::vec3(inverse_rot * glm::vec4(vec, 1)));
		rot = glm::rotate(rot,angle,new_vec);
	} else {
		rot = glm::rotate(rot,angle,vec);
	}

}
	
void  MovableGLM::MyScale(glm::vec3 scale)
{
	scl = glm::scale(scl,scale);
}

void MovableGLM::ZeroTrans()
{
	trans = glm::mat4(1);
	rot = glm::mat4(1);
	scl = glm::mat4(1);
}

glm::mat4 MovableGLM::get_trans()
{
	// return glm::mat4(trans);
	return trans;
}

glm::mat4 MovableGLM::get_rot()
{
	return glm::mat4(rot);
}

glm::mat4 MovableGLM::get_scl()
{
	return glm::mat4(scl);
}