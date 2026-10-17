"""Draw app icons with Qt; no external image assets are needed."""
from pathlib import Path
from PyQt5.QtGui import QImage, QPainter, QColor, QPen
from PyQt5.QtCore import QPointF

output = Path(__file__).resolve().parents[1] / 'web' / 'site'
output.mkdir(parents=True, exist_ok=True)
for size in (192, 512):
    image = QImage(size, size, QImage.Format_ARGB32)
    image.fill(QColor('#244263'))
    painter = QPainter(image)
    painter.setRenderHint(QPainter.Antialiasing)
    painter.scale(size / 192, size / 192)
    painter.setPen(QPen(QColor('#e9eff5'), 8))
    painter.drawLine(QPointF(48, 144), QPointF(85, 82))
    painter.drawLine(QPointF(85, 82), QPointF(123, 144))
    painter.drawLine(QPointF(40, 146), QPointF(139, 146))
    painter.setPen(QPen(QColor('#90b4d7'), 10))
    painter.drawLine(QPointF(63, 111), QPointF(117, 48))
    painter.setPen(QPen(QColor('#e9eff5'), 3))
    painter.drawLine(QPointF(117, 48), QPointF(148, 78))
    painter.setPen(QPen(QColor('#ff881d'), 1))
    painter.setBrush(QColor('#ff881d'))
    painter.drawEllipse(QPointF(148, 78), 10, 10)
    painter.end()
    assert image.save(str(output / f'icon-{size}.png'))
