#pragma once
#include <QAbstractTableModel>
#include <QTableView>

class BasicTable : public QAbstractTableModel {
public:
	BasicTable(const std::vector< std::vector<float> > & data , QObject *parent = 0) 
		: QAbstractTableModel(parent){
		Columns = data;
	}

	QVariant data(const QModelIndex& index, int role) const{
		if(role == Qt::DisplayRole)
		{
			return QString::number( Columns[index.row()][index.column()], 'g', 3 );
		}
		return QVariant::Invalid;
	}

	int rowCount(const QModelIndex& ) const{
		return (int)Columns.size();
	}

	int columnCount(const QModelIndex& ) const{
		return (int)Columns[0].size();
	}

	std::vector< std::vector<float> > Columns;
};

inline void showTable( std::vector< std::vector<float> > data, size_t count )
{
	static QTableView * tableview = new QTableView;

	data.resize(count);
	auto table = new BasicTable(data);

	tableview->setModel(table);
	tableview->setMinimumSize(400,900);
	tableview->show();
}

inline void showTable( const std::vector<float> & _data )
{
	std::vector< std::vector<float> > data;
	for(size_t i = 0; i < _data.size(); i++){
		std::vector<float> row;
		row.push_back(_data[i]);
		data.push_back(row);
	}

	static QTableView * tableview = new QTableView;
	auto table = new BasicTable(data);
	tableview->setModel(table);
	tableview->setMinimumSize(400,900);
	tableview->show();
}
