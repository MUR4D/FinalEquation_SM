double dmax4(double ay, double ay1, double ay2, double ay3)
{
    double w;
    w = ay;
    if (ay1 > w) w = ay1;
    if (ay2 > w) w = ay2;
    if (ay3 > w) w = ay3;
    return w;
}
