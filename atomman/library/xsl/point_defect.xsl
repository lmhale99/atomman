<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" 
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
  xmlns="http://www.w3.org/TR/xhtml1/strict">
  <xsl:output method="html" encoding="utf-8" indent="yes" />
  
  <xsl:template match="point-defect">
    <div>

      <style>
        .aslist {list-style-type: circle; list-style-position: inside; margin: 10px;}
      </style>

      <h1>Point defect parameter set</h1>
      <ul>
        <li><b><xsl:text>ID: </xsl:text></b><xsl:value-of select="id"/></li>
        <li><b><xsl:text>UUID4: </xsl:text></b><xsl:value-of select="key"/></li>
        <li><b><xsl:text>Family: </xsl:text></b>
          <xsl:choose>
            <xsl:when test="system-family-URL">
              <a href="{system-family-URL}"><xsl:value-of select="system-family"/></a>
            </xsl:when>
            <xsl:otherwise>
              <xsl:value-of select="system-family"/>
            </xsl:otherwise>
          </xsl:choose>
        </li>
        <xsl:for-each select="calculation-parameter">
          <li>
            <b><xsl:text>atomman.defect.point defect parameters:</xsl:text></b>
            <ul class="aslist">
              <li><xsl:text>ptd_type: "</xsl:text><xsl:value-of select="ptd_type"/><xsl:text>"</xsl:text></li>
              <xsl:if test="atype">
                <li><xsl:text>atype: "</xsl:text><xsl:value-of select="atype"/><xsl:text>"</xsl:text></li>
              </xsl:if>
              <li><xsl:text>pos: "</xsl:text><xsl:value-of select="pos"/><xsl:text>"</xsl:text></li>
              <xsl:if test="db_vect">
                <li><xsl:text>db_vect: "</xsl:text><xsl:value-of select="db_vect"/><xsl:text>"</xsl:text></li>
              </xsl:if>
              <li><xsl:text>scale: "</xsl:text><xsl:value-of select="scale"/><xsl:text>"</xsl:text></li>
            </ul>
          </li>
        </xsl:for-each>
      </ul>
    </div>

  </xsl:template>
</xsl:stylesheet>