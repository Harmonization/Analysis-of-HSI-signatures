# Анализ сигнатур HSI

Спектральной сигнатурой называется пиксель гиперспектрального изображения (HSI), представляющий собой вектор значений всех каналов HSI. Сигнатура может быть представлена в виде графика зависимости значений каналов от номера канала/длины волны. Такие графики имеют уникальную форму для разных материалов, и поэтому сигнатура гиперспектрального изображения может быть полезна в задачах классификации. 

Средним спектром (сигнатурой) называется вектор с усредненными значениями сигнатур HSI, лежащих в некоторой области (ROI). Такая сигнатура представляет собой усреднение всех материалов, лежащих в выбранной области. Если же все пиксели из данной области характеризуют один объект, то находится средняя сигнатура для данного объекта, форма графика которой будет более представительной, чем форма для случайно выбранного одиночного пикселя. 

Кроме сигнатур интерес также представляют спектральные индексы, представляющие собой одноканальное изображение, полученное математикой определенных каналов HSI. Одним из популярных классов индексов является (иногда нормализованная) разность двух выбранных каналов. Такие индексы могут эффективно выделять различные материалы (воду, дерево, почву, растения и тд), а также определять некоторые их свойства. Например вегетативные индексы (ВИ) позволяют эффективно сегментировать растительность на многоканальных изображениях. Эффективность сегментации определяется величиной разности каналов, которая чем больше по модулю, тем лучше выделяется материал на результирующем изображении индекса. Для исследования эффективности индексов полезно рассмотреть специальную матрицу, элементы которой определяются как: A[i,j] = MS[j] - MS[i], где MS - это средний спектр для некоторой области интереса (ROI). Из такой матрице можно получить каналы для индексов, наиболее подходящие для текущей задачи. 

# Описание программы

Данный код был написан для анализа сигнатур из изображений пшеницы, но может быть обобщен и для более широкого класса гиперспектральных изображений. Он представляет собой простой программный интерфейс для первичного анализа сигнатур и каналов HSI. Код интерфейса приведен в данном репозитории, а исполнительный файл программы лежит по ссылке: https://disk.yandex.ru/d/8AdPWpJrq67HAw. 

При открытии заголовочного файла с расширением .hdr откроется окно программы, в котором отобразится псевдо-RGB изображение, соответствующее открытому HSI. Нажатием правой кнопкой мыши на пиксель изображения в графике справа визуализируется спектральная сигнатура, соответствующая выбранному пикселю. Можно добавлять сигнатуры на график, очистить график или сохранить их в виде таблицы, нажав соответствующие кнопки. Также под изображением есть кнопка лупы, который позволяет приблизить выбранную область и выбрать требуемый пиксель наиболее точно. 

При нажатии кнопки "Выбрать ROI" окно программы меняется и появляются новые возможности. Растение на HSI отделяется с использованием ВИ, каналы которого задаются самим пользователем. После вычисления ВИ растение отделяется по порогу, задаваемому ползунком в интерфейсе. Меняя каналы индекса можно получить разные результаты - одни каналы оставляют только растение, другие оставляют только деревянные палки, третьи оставляют и то и другое. Не обязательно подбирать порог и каналы таким образом, чтобы удалить посторонние элементы со всего изображения. Достаточно выбрать один или несколько ROI и добиться удаления лишних материалов именно там. ROI выбираются по двум точкам с помощью правого клика мышью. В интерфейсе по умолчанию показывается все изображение целиком, но отметив соответствующий переключатель, можно визуализировать только ROI. 

С правой стороны интерфейса программы визуализированы график среднего спектра для выбранного ROI и соответствующая ему матрица разниц. По умолчанию они вычисляются для всего изображения, но с выбором ROI они будут автоматически пересчитаны для него. Чтобы вернуться к графикам предыдущего ROI (или всего изображения) нужно выбрать соответствующую вкладку. Средний спектр и матрица считаются только для части изображения, которая отделена спектральным индексом по выбранному порогу. Чтобы вычислить их для всего изображения, необходимо сделать порог отрицательным. Поскольку зачастую нужны не все каналы HSI, то есть возможность указать границы, между которыми будет визуализирована матрица разниц (текстовые поля под графиком). Значения среднего спектра и матрицы для всех отмеченных ROI можно сохранить в виде таблиц. 

![Иллюстрация к проекту](image/img_1.jpg)

