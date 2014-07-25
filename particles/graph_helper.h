#pragma once

#include "svgview.h"

template<typename GraphType>
inline QString toGraphvizFormat( GraphType & graph, QMap<int,QColor> nodeColors, QMap<int,QString> nodeLables ){
	QStringList out;
	out << "graph G{";
	out << "\t node [ fontcolor = black, color = white, style = filled ];";
	for(auto v : graph.vertices){
		QColor color = nodeColors[v];
		QString colorHex; colorHex.sprintf("#%02X%02X%02X", color.red(), color.green(), color.blue());
		out << ("\t" + QString("%1 [label = \"%2\", color = \"%3\"]").arg(v).arg(nodeLables[v]).arg(colorHex));
	}
	for(auto e : graph.GetEdgesSet()){
		QString elabel;
		auto e_weight = e.weight;
		if(e_weight != 1) elabel = QString::number(e_weight);

		out << ("\t\"" + QString::number(e.index) + "\" -- \"" + QString::number(e.target) + QString("\" [ label=\"%1\" ]").arg(elabel));
	}
	out << "fontsize = 20;";
	out << "}";
	return out.join("\n");
}

#include <QProcess>
inline QString buildSVG( QString graphVizText )
{
	// Locate 'dot' application
	static QString dotPath;
	if(!dotPath.size())
	{
		QString path = std::getenv("PATH");
		QStringList path_folders = path.split(";");
		path_folders.sort();
		
		for(auto folder : path_folders){
			QString exe = folder.replace("\\", "/") + "/dot.exe";
			QFileInfo info( exe );
			if(info.exists()){
				dotPath = exe;
			}
		}
	}

	auto p = new QProcess;        
	p->start( dotPath, QStringList() << "-Tsvg" );
	p->write( qPrintable(graphVizText) );
	p->closeWriteChannel();
	p->waitForFinished();

	QString svgoutput (p->readAllStandardOutput());
	return svgoutput;
}
