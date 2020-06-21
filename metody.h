#ifndef METODY_H
#define METODY_H

#include <QMainWindow>

namespace Ui {
class metody;
}

class metody : public QMainWindow
{
    Q_OBJECT

public:
    explicit metody(QWidget *parent = nullptr);
    ~metody();

private slots:

    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

private:
    Ui::metody *ui;
};

#endif // METODY_H
