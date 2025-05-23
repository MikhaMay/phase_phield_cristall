### Proof of concept

- Выбрана цель — визуализация 3d или добавление в модель упругости
- Программа считает 1d на уровне старой реализации, но на порядки быстрее
- Программа корректно считает 3d или запрограммированы базовые формулы для упругости (из статьи)

### Особенности реализации

#### Способы ускорения:

- ~~параллелить на cpu~~
- ~~попробовать на gpu (пойдет в защиту проекта по параллельным вычислениям), возможно на python~~
  параллелить не вариант — сетка мала, накладных расходов на потоки больше
- поднять виртуальную машину с большим числом cpu / с gpu (???)


#### Задачи

[x] Настроить vs code
- создать рабочую область, настроить под C++
- подключить систему контроля версий
- первый удачный запуск старого кода
- запушить на гитхаб

[x] Оптимизация
- ~~распараллелить (пока на cpu)~~  
    omp не дало результатов — слишком малая область сетки; зато помогла команда компиляции, включающая оптимизацию + автовекторизация

[x] Рефакторинг
- разбить реализацию с большого файла на модули
- задавать все параметры через конфиг ~~(например, yaml)~~, а не в коде (включая тип ГУ из пресета)  
    в качестве конфига выбрал .h файл, хотелось уметь писать формулы в параметрах и рассчитывать одни параметры из других
- удобное сохранение бинарных файлов с результатами расчёта  
    результаты запусков с отличающимися параметрами должны лежать отдельно и не должны перетирать друг друга

[x] Проверить все формулы на корректность (сравнить с текстом курсовой)

[x] Потестировать код локально на больших объемах, вывезет ли макбук

[x] Попробовать сделать 2d версию

[x] ~~Потестировать на gpu~~

[x] Решить поднимать ли виртуалку со слабой gpu или с большим количеством cpu, бюджет грантовые 20к рублей на яндекс клауд, либо поднять в сторонних сервисах

[x] Посчиать 2d плавление в половине области

[ ] Посчиать 2d кристализацию в половине области

[x] Подправить препринт

[x] Подправить теорию по sHPFC модели

[ ] Написать схему для sHPFC

[ ] Добавить в курсач раздел про sHPFC кристалл

[ ] Найти статью, запрограммировать формулу для упругости
