#include <QPixmap>
#include <stack>
#include "LFDWidget.h"
#include "ui_LFDWidget.h"

#include "icosahedron.h"

LFDWidget::LFDWidget(SurfaceMesh::SurfaceMeshModel *model, QWidget *parent) : QWidget(parent), ui(new Ui::LFDWidget)
{
    ui->setupUi(this);
    ui->mainlayout->addWidget(new MyGLWidget(model, parent));
}

LFDWidget::~LFDWidget()
{
    delete ui;
}

void MyGLWidget::initializeGL()
{
	this->setFixedSize(QSize(1280, 720));

	initializeOpenGLFunctions();
	glClearColor(0, 0, 0, 1);

	glDisable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);

	// Setup lights and material
	GLfloat ambientLightColor[] = { 0.2f, 0.2f, 0.2f, 1 };
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLightColor);

	GLfloat diffuseLightColor[] = { 0.9f, 0.9f, 0.9f, 1 };
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLightColor);

	GLfloat specularLightColor[] = { 0.95f, 0.95f, 0.95f, 1 };
	glLightfv(GL_LIGHT0, GL_SPECULAR, specularLightColor);

	float posLight0[] = { 3, 3, 3, 0 };
	glLightfv(GL_LIGHT0, GL_POSITION, posLight0);

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	// Specular lighting
	float specReflection[] = { 0.8f, 0.8f, 0.8f, 1.0f };
	glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
	glMateriali(GL_FRONT, GL_SHININESS, 56);

	this->isReady = true;
}

void MyGLWidget::paintGL()
{
	if (!this->isReady) return;
	if (this->isDone) return;

	// Sample uniformly on sphere:
	int sample_level = 4;
	auto samples = icosahedron::sample(sample_level);

	// Resolution
	int width = 128;
	int cols = std::floor(this->width() / width);
	int rows = std::floor(this->height() / width);

	// Load geometry
	static GLuint VertexVBOID = -1;
	if (VertexVBOID == -1)
	{
		// Update
		{
			model->updateBoundingBox();
			model->update_face_normals();
			model->update_vertex_normals();
		}

		// Normalize
		{
			Eigen::AlignedBox3d bbox = model->bbox();
			Eigen::Vector3d offset = bbox.center();
			Eigen::Vector3d s = bbox.diagonal();
			double scale = std::max(s.x(), std::max(s.y(), s.z()));
			SurfaceMesh::Vector3VertexProperty points = model->vertex_coordinates();
			for(auto v : model->vertices()){
				auto & p = points[v];
				p = (p - offset) / scale;
			}
		}

		std::vector<GLVertexd> all_vertices;
		for (auto v : model->vertices()) {
			auto p = model->vertex_coordinates()[v];
			auto n = model->vertex_normals()[v];
			all_vertices.push_back(GLVertexd(p.x(), p.y(), p.z(), n.x(), n.y(), n.z()));
		}

		std::vector<GLVertexd> vertices;
		for (auto f : model->faces()){
			for (auto v : model->vertices(f))
				vertices.push_back(all_vertices[v.idx()]);
		}

		// Vertices
		glGenBuffers(1, &VertexVBOID);
		glBindBuffer(GL_ARRAY_BUFFER, VertexVBOID);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLVertexd) * vertices.size(), &vertices[0].x, GL_STATIC_DRAW);
	}

	auto drawView = [&](int x, int y, int size, RenderOption & render){
		// Setup viewport
		glViewport(x, y, size, size);
		glScissor(x, y, size, size);

		// Background
		if (render.isBackground)
		{
			// Setup 2D
			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
			glLoadIdentity();
			glOrtho(0, size, size, 0, 0.0, -1.0);
			glMatrixMode(GL_MODELVIEW);
			glPushMatrix();
			glLoadIdentity();

			// Draw gradient quad
			glDisable(GL_LIGHTING);
			glBegin(GL_QUADS);
			glColor3d(0, 0, 0); glVertex2d(0, 0); glVertex2d(size, 0);
			glColor3d(0.15, 0.15, 0.15); glVertex2d(size, size); glVertex2d(0, size);
			glEnd();
			glEnable(GL_LIGHTING);

			// End 2D
			glMatrixMode(GL_PROJECTION);
			glPopMatrix();
			glMatrixMode(GL_MODELVIEW);
			glPopMatrix();

			// Background has no depth
			glClear(GL_DEPTH_BUFFER_BIT);
		}

		// Camera
		qglviewer::Camera cam;
		cam.setType(qglviewer::Camera::ORTHOGRAPHIC);
		cam.setScreenWidthAndHeight(size, size);
		cam.setSceneRadius(20.0f);
		cam.setUpVector(qglviewer::Vec(0, 0, 1));
		cam.setSceneCenter(qglviewer::Vec(0, 0, 0));
		cam.setPosition(render.cameraPos);
		cam.setViewDirection(-render.cameraPos);
		cam.loadProjectionMatrix();
		cam.loadModelViewMatrix();

		// Color
		if(render.isColor) glColor3f(0, 1, 0);
		else glColor3f(1, 1, 1);

		// Lights
		if (render.isLights) glEnable(GL_LIGHTING);
		else glDisable(GL_LIGHTING);

		// Draw
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		glBindBuffer(GL_ARRAY_BUFFER, VertexVBOID);
		glVertexPointer(3, GL_DOUBLE, sizeof(GLVertexd), (void*)offsetof(GLVertexd, x));
		glNormalPointer(GL_DOUBLE, sizeof(GLVertexd), (void*)offsetof(GLVertexd, nx));
		glDrawArrays(GL_TRIANGLES, 0, model->n_faces() * 3);

		// Visualize cameras
		if (render.isShowCameras)
		{
			glDisable(GL_LIGHTING);
			glPointSize(3);
			glColor3d(1, 1, 1);
			glBegin(GL_POINTS);
			for (auto p : samples) { glColor3d(p.x, p.y, p.z); glVertex3d(p.x, p.y, p.z); }
			glEnd();
		}
	};

	// Create render jobs
	std::stack < RenderOption > views;
	for (auto point : samples){
		double scaling = 1.8;
		auto cameraPos = scaling * qglviewer::Vec(point.x, point.y, point.z);
		RenderOption render_option(cameraPos);

		render_option.isBackground = render_option.isShowCameras = render_option.isColor = false;
		render_option.isLights = false;

		views.push(render_option);
	}

	// Convert samples into batches to draw per frame
	std::vector < std::vector < RenderOption > > batch(1);
	int batch_capacity = cols * rows;
	while (!views.empty()){
		auto v = views.top();
		views.pop();
		if (batch.back().size() >= batch_capacity) batch.push_back(std::vector < RenderOption >());
		batch.back().push_back(v);
	}

	QVector < QImage > silhouettes;

	for (auto & b : batch)
	{
		int w = this->width(), h = this->height();

		glViewport(0, 0, w, h);
		glScissor(0, 0, w, h);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		int idx = 0;
		for (int y = 0; y < rows; y++){
			for (int x = 0; x < cols; x++){
				int dx = x * width, dy = y * width;
				drawView(dx, dy, width, b[idx++]);

				// Break on partial batches;
				if (idx == b.size()){
					x = cols;
					y = rows;
				}
			}
		}

		// Grab buffer
		QImage res = QImage(w, h, QImage::Format_ARGB32);
		glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, res.bits());
		if(b.front().isColor) res = res.rgbSwapped();

		// Split into images
		auto splitAll = [&]( const QImage & img, int count ){
			QVector<QImage> pieces;
			int idx = 0;
			for (int y = 0; y < rows; y++){
				for (int x = 0; x < cols; x++){
					int dx = x * width, dy = y * width;
					pieces.push_back(img.copy( QRect(dx,dy,width,width) ));
					idx++;
					if (idx == count){
						x = cols;
						y = rows;
					}
				}
			}
			return pieces;
		};

		silhouettes += splitAll(res, b.size());
	}

	// Save to files
	for (int i = 0; i < silhouettes.size(); i++){
		QString number; number.sprintf("%05d", i);
		silhouettes[i].save( QString("s%1.png").arg(number) );
	}

	isDone = true;
}
